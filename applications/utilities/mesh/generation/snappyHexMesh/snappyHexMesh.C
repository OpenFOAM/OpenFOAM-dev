/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    snappyHexMesh

Description
    Automatic split hex mesher. Refines and snaps to surface.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "snappyRefineDriver.H"
#include "snappySnapDriver.H"
#include "snappyLayerDriver.H"
#include "searchableSurfaces.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "refinementRegions.H"
#include "decompositionMethod.H"
#include "noDecomp.H"
#include "fvMeshDistribute.H"
#include "wallPolyPatch.H"
#include "refinementParameters.H"
#include "snapParameters.H"
#include "layerParameters.H"
#include "faceSet.H"
#include "motionSmoother.H"
#include "polyTopoChange.H"
#include "cellModeller.H"
#include "uindirectPrimitivePatch.H"
#include "surfZoneIdentifierList.H"
#include "UnsortedMeshedSurface.H"
#include "MeshedSurface.H"
#include "globalIndex.H"
#include "IOmanip.H"
#include "fvMeshTools.H"
#include "systemDict.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Convert size (as fraction of defaultCellSize) to refinement level
label sizeCoeffToRefinement
(
    const scalar level0Coeff,   // ratio of hex cell size v.s. defaultCellSize
    const scalar sizeCoeff
)
{
     return round(::log(level0Coeff/sizeCoeff)/::log(2));
}


autoPtr<refinementSurfaces> createRefinementSurfaces
(
    const searchableSurfaces& allGeometry,
    const dictionary& surfacesDict,
    const dictionary& shapeControlDict,
    const label gapLevelIncrement,
    const scalar level0Coeff
)
{
    autoPtr<refinementSurfaces> surfacePtr;

    // Count number of surfaces.
    label surfI = 0;
    forAll(allGeometry.names(), geomI)
    {
        const word& geomName = allGeometry.names()[geomI];

        if (surfacesDict.found(geomName))
        {
            surfI++;
        }
    }

    labelList surfaces(surfI);
    wordList names(surfI);
    PtrList<surfaceZonesInfo> surfZones(surfI);

    labelList regionOffset(surfI);

    labelList globalMinLevel(surfI, 0);
    labelList globalMaxLevel(surfI, 0);
    labelList globalLevelIncr(surfI, 0);
    PtrList<dictionary> globalPatchInfo(surfI);
    List<Map<label>> regionMinLevel(surfI);
    List<Map<label>> regionMaxLevel(surfI);
    List<Map<label>> regionLevelIncr(surfI);
    List<Map<scalar>> regionAngle(surfI);
    List<Map<autoPtr<dictionary>>> regionPatchInfo(surfI);

    HashSet<word> unmatchedKeys(surfacesDict.toc());

    surfI = 0;
    forAll(allGeometry.names(), geomI)
    {
        const word& geomName = allGeometry.names()[geomI];

        const entry* ePtr = surfacesDict.lookupEntryPtr(geomName, false, true);

        if (ePtr)
        {
            const dictionary& shapeDict = ePtr->dict();
            unmatchedKeys.erase(ePtr->keyword());

            names[surfI] = geomName;
            surfaces[surfI] = geomI;

            const searchableSurface& surface = allGeometry[geomI];

            // Find the index in shapeControlDict
            // Invert surfaceCellSize to get the refinementLevel

            const word scsFuncName =
                shapeDict.lookup("surfaceCellSizeFunction");
            const dictionary& scsDict =
                shapeDict.optionalSubDict(scsFuncName + "Coeffs");

            const scalar surfaceCellSize =
                scsDict.lookup<scalar>("surfaceCellSizeCoeff");

            const label refLevel = sizeCoeffToRefinement
            (
                level0Coeff,
                surfaceCellSize
            );

            globalMinLevel[surfI] = refLevel;
            globalMaxLevel[surfI] = refLevel;
            globalLevelIncr[surfI] = gapLevelIncrement;

            // Surface zones
            surfZones.set(surfI, new surfaceZonesInfo(surface, shapeDict));


            // Global perpendicular angle
            if (shapeDict.found("patchInfo"))
            {
                globalPatchInfo.set
                (
                    surfI,
                    shapeDict.subDict("patchInfo").clone()
                );
            }


            // Per region override of patchInfo

            if (shapeDict.found("regions"))
            {
                const dictionary& regionsDict = shapeDict.subDict("regions");
                const wordList& regionNames =
                    allGeometry[surfaces[surfI]].regions();

                forAll(regionNames, regionI)
                {
                    if (regionsDict.found(regionNames[regionI]))
                    {
                        // Get the dictionary for region
                        const dictionary& regionDict = regionsDict.subDict
                        (
                            regionNames[regionI]
                        );

                        if (regionDict.found("patchInfo"))
                        {
                            regionPatchInfo[surfI].insert
                            (
                                regionI,
                                regionDict.subDict("patchInfo").clone()
                            );
                        }
                    }
                }
            }

            // Per region override of cellSize
            if (shapeDict.found("regions"))
            {
                const dictionary& shapeControlRegionsDict =
                    shapeDict.subDict("regions");
                const wordList& regionNames =
                    allGeometry[surfaces[surfI]].regions();

                forAll(regionNames, regionI)
                {
                    if (shapeControlRegionsDict.found(regionNames[regionI]))
                    {
                        const dictionary& shapeControlRegionDict =
                            shapeControlRegionsDict.subDict
                            (
                                regionNames[regionI]
                            );

                        const word scsFuncName =
                            shapeControlRegionDict.lookup
                            (
                                "surfaceCellSizeFunction"
                            );
                        const dictionary& scsDict =
                            shapeControlRegionDict.subDict
                            (
                                scsFuncName + "Coeffs"
                            );

                        const scalar surfaceCellSize =
                                scsDict.lookup<scalar>("surfaceCellSizeCoeff");

                        const label refLevel = sizeCoeffToRefinement
                        (
                            level0Coeff,
                            surfaceCellSize
                        );

                        regionMinLevel[surfI].insert(regionI, refLevel);
                        regionMaxLevel[surfI].insert(regionI, refLevel);
                        regionLevelIncr[surfI].insert(regionI, 0);
                    }
                }
            }

            surfI++;
        }
    }

    // Calculate local to global region offset
    label nRegions = 0;

    forAll(surfaces, surfI)
    {
        regionOffset[surfI] = nRegions;
        nRegions += allGeometry[surfaces[surfI]].regions().size();
    }

    // Rework surface specific information into information per global region
    labelList minLevel(nRegions, 0);
    labelList maxLevel(nRegions, 0);
    labelList gapLevel(nRegions, -1);
    PtrList<dictionary> patchInfo(nRegions);

    forAll(globalMinLevel, surfI)
    {
        label nRegions = allGeometry[surfaces[surfI]].regions().size();

        // Initialise to global (i.e. per surface)
        for (label i = 0; i < nRegions; i++)
        {
            label globalRegionI = regionOffset[surfI] + i;
            minLevel[globalRegionI] = globalMinLevel[surfI];
            maxLevel[globalRegionI] = globalMaxLevel[surfI];
            gapLevel[globalRegionI] =
                maxLevel[globalRegionI]
              + globalLevelIncr[surfI];

            if (globalPatchInfo.set(surfI))
            {
                patchInfo.set
                (
                    globalRegionI,
                    globalPatchInfo[surfI].clone()
                );
            }
        }

        // Overwrite with region specific information
        forAllConstIter(Map<label>, regionMinLevel[surfI], iter)
        {
            label globalRegionI = regionOffset[surfI] + iter.key();

            minLevel[globalRegionI] = iter();
            maxLevel[globalRegionI] = regionMaxLevel[surfI][iter.key()];
            gapLevel[globalRegionI] =
                maxLevel[globalRegionI]
              + regionLevelIncr[surfI][iter.key()];
        }

        const Map<autoPtr<dictionary>>& localInfo = regionPatchInfo[surfI];
        forAllConstIter(Map<autoPtr<dictionary>>, localInfo, iter)
        {
            label globalRegionI = regionOffset[surfI] + iter.key();
            patchInfo.set(globalRegionI, iter()().clone());
        }
    }

    surfacePtr.set
    (
        new refinementSurfaces
        (
            allGeometry,
            surfaces,
            names,
            surfZones,
            regionOffset,
            minLevel,
            maxLevel,
            gapLevel,
            scalarField(nRegions, -great),  // perpendicularAngle,
            patchInfo
        )
    );


    const refinementSurfaces& rf = surfacePtr();

    // Determine maximum region name length
    label maxLen = 0;
    forAll(rf.surfaces(), surfI)
    {
        label geomI = rf.surfaces()[surfI];
        const wordList& regionNames = allGeometry.regionNames()[geomI];
        forAll(regionNames, regionI)
        {
            maxLen = Foam::max(maxLen, label(regionNames[regionI].size()));
        }
    }


    Info<< setw(maxLen) << "Region"
        << setw(10) << "Min Level"
        << setw(10) << "Max Level"
        << setw(10) << "Gap Level" << nl
        << setw(maxLen) << "------"
        << setw(10) << "---------"
        << setw(10) << "---------"
        << setw(10) << "---------" << endl;

    forAll(rf.surfaces(), surfI)
    {
        label geomI = rf.surfaces()[surfI];

        Info<< rf.names()[surfI] << ':' << nl;

        const wordList& regionNames = allGeometry.regionNames()[geomI];

        forAll(regionNames, regionI)
        {
            label globalI = rf.globalRegion(surfI, regionI);

            Info<< setw(maxLen) << regionNames[regionI]
                << setw(10) << rf.minLevel()[globalI]
                << setw(10) << rf.maxLevel()[globalI]
                << setw(10) << rf.gapLevel()[globalI] << endl;
        }
    }


    return surfacePtr;
}


void extractSurface
(
    const polyMesh& mesh,
    const Time& runTime,
    const labelHashSet& includePatches,
    const fileName& outFileName
)
{
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    // Collect sizes. Hash on names to handle local-only patches (e.g.
    //  processor patches)
    HashTable<label> patchSize(1000);
    label nFaces = 0;
    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const polyPatch& pp = bMesh[iter.key()];
        patchSize.insert(pp.name(), pp.size());
        nFaces += pp.size();
    }
    Pstream::mapCombineGather(patchSize, plusEqOp<label>());


    // Allocate zone/patch for all patches
    HashTable<label> compactZoneID(1000);
    forAllConstIter(HashTable<label>, patchSize, iter)
    {
        label sz = compactZoneID.size();
        compactZoneID.insert(iter.key(), sz);
    }
    Pstream::mapCombineScatter(compactZoneID);


    // Rework HashTable into labelList just for speed of conversion
    labelList patchToCompactZone(bMesh.size(), -1);
    forAllConstIter(HashTable<label>, compactZoneID, iter)
    {
        label patchi = bMesh.findPatchID(iter.key());
        if (patchi != -1)
        {
            patchToCompactZone[patchi] = iter();
        }
    }

    // Collect faces on zones
    DynamicList<label> faceLabels(nFaces);
    DynamicList<label> compactZones(nFaces);
    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const polyPatch& pp = bMesh[iter.key()];
        forAll(pp, i)
        {
            faceLabels.append(pp.start()+i);
            compactZones.append(patchToCompactZone[pp.index()]);
        }
    }

    // Addressing engine for all faces
    uindirectPrimitivePatch allBoundary
    (
        UIndirectList<face>(mesh.faces(), faceLabels),
        mesh.points()
    );


    // Find correspondence to master points
    labelList pointToGlobal;
    labelList uniqueMeshPoints;
    autoPtr<globalIndex> globalNumbers = mesh.globalData().mergePoints
    (
        allBoundary.meshPoints(),
        allBoundary.meshPointMap(),
        pointToGlobal,
        uniqueMeshPoints
    );

    // Gather all unique points on master
    List<pointField> gatheredPoints(Pstream::nProcs());
    gatheredPoints[Pstream::myProcNo()] = pointField
    (
        mesh.points(),
        uniqueMeshPoints
    );
    Pstream::gatherList(gatheredPoints);

    // Gather all faces
    List<faceList> gatheredFaces(Pstream::nProcs());
    gatheredFaces[Pstream::myProcNo()] = allBoundary.localFaces();
    forAll(gatheredFaces[Pstream::myProcNo()], i)
    {
        inplaceRenumber(pointToGlobal, gatheredFaces[Pstream::myProcNo()][i]);
    }
    Pstream::gatherList(gatheredFaces);

    // Gather all ZoneIDs
    List<labelList> gatheredZones(Pstream::nProcs());
    gatheredZones[Pstream::myProcNo()] = move(compactZones);
    Pstream::gatherList(gatheredZones);

    // On master combine all points, faces, zones
    if (Pstream::master())
    {
        pointField allPoints = ListListOps::combine<pointField>
        (
            gatheredPoints,
            accessOp<pointField>()
        );
        gatheredPoints.clear();

        faceList allFaces = ListListOps::combine<faceList>
        (
            gatheredFaces,
            accessOp<faceList>()
        );
        gatheredFaces.clear();

        labelList allZones = ListListOps::combine<labelList>
        (
            gatheredZones,
            accessOp<labelList>()
        );
        gatheredZones.clear();


        // Zones
        surfZoneIdentifierList surfZones(compactZoneID.size());
        forAllConstIter(HashTable<label>, compactZoneID, iter)
        {
            surfZones[iter()] = surfZoneIdentifier(iter.key(), iter());
            Info<< "surfZone " << iter()  <<  " : " << surfZones[iter()].name()
                << endl;
        }

        UnsortedMeshedSurface<face> unsortedFace
        (
            move(allPoints),
            move(allFaces),
            move(allZones),
            move(surfZones)
        );


        MeshedSurface<face> sortedFace(unsortedFace);

        Info<< "Writing merged surface to "
            << runTime.globalPath()/outFileName << endl;

        sortedFace.write(runTime.globalPath()/outFileName);
    }
}


// Check writing tolerance before doing any serious work
scalar getMergeDistance(const polyMesh& mesh, const scalar mergeTol)
{
    const boundBox& meshBb = mesh.bounds();
    scalar mergeDist = mergeTol * meshBb.mag();

    Info<< nl
        << "Overall mesh bounding box  : " << meshBb << nl
        << "Relative tolerance         : " << mergeTol << nl
        << "Absolute matching distance : " << mergeDist << nl
        << endl;

    // check writing tolerance
    if (mesh.time().writeFormat() == IOstream::ASCII)
    {
        const scalar writeTol = std::pow
        (
            scalar(10.0),
            -scalar(IOstream::defaultPrecision())
        );

        if (mergeTol < writeTol)
        {
            FatalErrorInFunction
                << "Your current settings specify ASCII writing with "
                << IOstream::defaultPrecision() << " digits precision." << nl
                << "Your merging tolerance (" << mergeTol
                << ") is finer than this." << nl
                << "Change to binary writeFormat, "
                << "or increase the writePrecision" << endl
                << "or adjust the merge tolerance (mergeTol)."
                << exit(FatalError);
        }
    }

    return mergeDist;
}


void removeZeroSizedPatches(fvMesh& mesh)
{
    // Remove any zero-sized ones. Assumes
    // - processor patches are already only there if needed
    // - all other patches are available on all processors
    // - but coupled ones might still be needed, even if zero-size
    //   (e.g. processorCyclic)
    // See also logic in createPatch.
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList oldToNew(pbm.size(), -1);
    label newPatchi = 0;
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (!isA<processorPolyPatch>(pp))
        {
            if
            (
                isA<coupledPolyPatch>(pp)
             || returnReduce(pp.size(), sumOp<label>())
            )
            {
                // Coupled (and unknown size) or uncoupled and used
                oldToNew[patchi] = newPatchi++;
            }
        }
    }

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            oldToNew[patchi] = newPatchi++;
        }
    }


    const label nKeepPatches = newPatchi;

    // Shuffle unused ones to end
    if (nKeepPatches != pbm.size())
    {
        Info<< endl
            << "Removing zero-sized patches:" << endl << incrIndent;

        forAll(oldToNew, patchi)
        {
            if (oldToNew[patchi] == -1)
            {
                Info<< indent << pbm[patchi].name()
                    << " type " << pbm[patchi].type()
                    << " at position " << patchi << endl;
                oldToNew[patchi] = newPatchi++;
            }
        }
        Info<< decrIndent;

        fvMeshTools::reorderPatches(mesh, oldToNew, nKeepPatches, true);
        Info<< endl;
    }
}


// Write mesh and additional information
void writeMesh
(
    const string& msg,
    const meshRefinement& meshRefiner,
    const meshRefinement::debugType debugLevel,
    const meshRefinement::writeType writeLevel
)
{
    const fvMesh& mesh = meshRefiner.mesh();

    meshRefiner.printMeshInfo(debugLevel, msg);
    Info<< "Writing mesh to time " << meshRefiner.timeName() << endl;

    meshRefiner.write
    (
        debugLevel,
        meshRefinement::writeType(writeLevel | meshRefinement::WRITEMESH),
        mesh.time().path()/meshRefiner.timeName()
    );
    Info<< "Wrote mesh in = "
        << mesh.time().cpuTimeIncrement() << " s." << endl;
}


int main(int argc, char *argv[])
{
    #include "addOverwriteOption.H"
    Foam::argList::addBoolOption
    (
        "checkGeometry",
        "check all surface geometry for quality"
    );
    Foam::argList::addOption
    (
        "surfaceSimplify",
        "boundBox",
        "simplify the surface using snappyHexMesh starting from a boundBox"
    );
    Foam::argList::addOption
    (
        "patches",
        "(patch0 .. patchN)",
        "only triangulate selected patches (wildcards supported)"
    );
    Foam::argList::addOption
    (
        "outFile",
        "file",
        "name of the file to save the simplified surface to"
    );
    #include "addDictOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();

    const bool overwrite = args.optionFound("overwrite");
    const bool checkGeometry = args.optionFound("checkGeometry");
    const bool surfaceSimplify = args.optionFound("surfaceSimplify");

    autoPtr<fvMesh> meshPtr;

    {
        Foam::Info
            << "Create mesh for time = "
            << runTime.timeName() << Foam::nl << Foam::endl;

        meshPtr.set
        (
            new fvMesh
            (
                Foam::IOobject
                (
                    Foam::fvMesh::defaultRegion,
                    runTime.timeName(),
                    runTime,
                    Foam::IOobject::MUST_READ
                )
            )
        );
    }

    fvMesh& mesh = meshPtr();

    Info<< "Read mesh in = "
        << runTime.cpuTimeIncrement() << " s" << endl;

    // Check patches and faceZones are synchronised
    mesh.boundaryMesh().checkParallelSync(true);
    meshRefinement::checkCoupledFaceZones(mesh);


    // Read meshing dictionary
    const IOdictionary meshDict(systemDict("snappyHexMeshDict", args, mesh));


    // all surface geometry
    const dictionary& geometryDict = meshDict.subDict("geometry");

    // refinement parameters
    const dictionary& refineDict = meshDict.subDict("castellatedMeshControls");

    // mesh motion and mesh quality parameters
    const dictionary& motionDict = meshDict.subDict("meshQualityControls");

    // snap-to-surface parameters
    const dictionary& snapDict = meshDict.subDict("snapControls");

    // layer addition parameters
    const dictionary& layerDict = meshDict.subDict("addLayersControls");

    // absolute merge distance
    const scalar mergeDist = getMergeDistance
    (
        mesh,
        meshDict.lookup<scalar>("mergeTolerance")
    );

    const Switch keepPatches(meshDict.lookupOrDefault("keepPatches", false));



    // Read decomposePar dictionary
    dictionary decomposeDict;
    {
        if (Pstream::parRun())
        {
            decomposeDict = IOdictionary
            (
                IOobject
                (
                    "decomposeParDict",
                    runTime.system(),
                    mesh,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            );
        }
        else
        {
            decomposeDict.add("method", "none");
            decomposeDict.add("numberOfSubdomains", 1);
        }
    }


    // Debug
    // ~~~~~

    // Set debug level
    meshRefinement::debugType debugLevel = meshRefinement::debugType
    (
        meshDict.lookupOrDefault<label>
        (
            "debug",
            0
        )
    );
    {
        wordList flags;
        if (meshDict.readIfPresent("debugFlags", flags))
        {
            debugLevel = meshRefinement::debugType
            (
                meshRefinement::readFlags
                (
                    meshRefinement::IOdebugTypeNames,
                    flags
                )
            );
        }
    }
    if (debugLevel > 0)
    {
        meshRefinement::debug   = debugLevel;
        snappyRefineDriver::debug = debugLevel;
        snappySnapDriver::debug   = debugLevel;
        snappyLayerDriver::debug  = debugLevel;
    }

    // Set file writing level
    {
        wordList flags;
        if (meshDict.readIfPresent("writeFlags", flags))
        {
            meshRefinement::writeLevel
            (
                meshRefinement::writeType
                (
                    meshRefinement::readFlags
                    (
                        meshRefinement::IOwriteTypeNames,
                        flags
                    )
                )
            );
        }
    }

    // Set output level
    {
        wordList flags;
        if (meshDict.readIfPresent("outputFlags", flags))
        {
            meshRefinement::outputLevel
            (
                meshRefinement::outputType
                (
                    meshRefinement::readFlags
                    (
                        meshRefinement::IOoutputTypeNames,
                        flags
                    )
                )
            );
        }
    }


    // Read geometry
    // ~~~~~~~~~~~~~

    searchableSurfaces allGeometry
    (
        IOobject
        (
            "abc",
            mesh.time().constant(),
            searchableSurface::geometryDir(mesh.time()),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        geometryDict,
        meshDict.lookupOrDefault("singleRegionName", true)
    );


    // Read refinement surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    autoPtr<refinementSurfaces> surfacesPtr;

    Info<< "Reading refinement surfaces..." << endl;

    if (surfaceSimplify)
    {
        IOdictionary foamyHexMeshDict
        (
           IOobject
           (
                "foamyHexMeshDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
           )
        );

        const dictionary& conformationDict =
            foamyHexMeshDict.subDict("surfaceConformation").subDict
            (
                "geometryToConformTo"
            );

        const dictionary& motionDict =
            foamyHexMeshDict.subDict("motionControl");

        const dictionary& shapeControlDict =
            motionDict.subDict("shapeControlFunctions");

        // Calculate current ratio of hex cells v.s. wanted cell size
        const scalar defaultCellSize =
            motionDict.lookup<scalar>("defaultCellSize");

        const scalar initialCellSize = ::pow(meshPtr().V()[0], 1.0/3.0);

        // Info<< "Wanted cell size  = " << defaultCellSize << endl;
        // Info<< "Current cell size = " << initialCellSize << endl;
        // Info<< "Fraction          = " << initialCellSize/defaultCellSize
        //    << endl;

        surfacesPtr =
            createRefinementSurfaces
            (
                allGeometry,
                conformationDict,
                shapeControlDict,
                refineDict.lookupOrDefault("gapLevelIncrement", 0),
                initialCellSize/defaultCellSize
            );
    }
    else
    {
        surfacesPtr.set
        (
            new refinementSurfaces
            (
                allGeometry,
                refineDict.found("refinementSurfaces")
              ? refineDict.subDict("refinementSurfaces")
              : dictionary::null,
                refineDict.lookupOrDefault("gapLevelIncrement", 0)
            )
        );

        Info<< "Read refinement surfaces in = "
            << mesh.time().cpuTimeIncrement() << " s" << nl << endl;
    }

    refinementSurfaces& surfaces = surfacesPtr();


    // Checking only?

    if (checkGeometry)
    {
        // Extract patchInfo
        List<wordList> patchTypes(allGeometry.size());

        const PtrList<dictionary>& patchInfo = surfaces.patchInfo();
        const labelList& surfaceGeometry = surfaces.surfaces();
        forAll(surfaceGeometry, surfI)
        {
            label geomI = surfaceGeometry[surfI];
            const wordList& regNames = allGeometry.regionNames()[geomI];

            patchTypes[geomI].setSize(regNames.size());
            forAll(regNames, regionI)
            {
                label globalRegionI = surfaces.globalRegion(surfI, regionI);

                if (patchInfo.set(globalRegionI))
                {
                    patchTypes[geomI][regionI] =
                        word(patchInfo[globalRegionI].lookup("type"));
                }
                else
                {
                    patchTypes[geomI][regionI] = wallPolyPatch::typeName;
                }
            }
        }

        // Write some stats
        allGeometry.writeStats(patchTypes, Info);
        // Check topology
        allGeometry.checkTopology(true);
        // Check geometry
        allGeometry.checkGeometry
        (
            100.0,      // max size ratio
            1e-9,       // intersection tolerance
            0.01,       // min triangle quality
            true
        );

        return 0;
    }



    // Read refinement shells
    // ~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Reading refinement regions..." << endl;
    refinementRegions shells
    (
        allGeometry,
        refineDict.found("refinementRegions")
      ? refineDict.subDict("refinementRegions")
      : dictionary::null
    );
    Info<< "Read refinement regions in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;


    Info<< "Setting refinement level of surface to be consistent"
        << " with shells.." << endl;
    surfaces.setMinLevelFields(shells);
    Info<< "Checked shell refinement in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;


    // Read feature meshes
    // ~~~~~~~~~~~~~~~~~~~

    Info<< "Reading features..." << endl;
    refinementFeatures features
    (
        mesh,
        refineDict.found("features")
      ? refineDict.lookup("features")
      : PtrList<dictionary>()
    );
    Info<< "Read features in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;



    // Refinement engine
    // ~~~~~~~~~~~~~~~~~

    Info<< nl
        << "Determining initial surface intersections" << nl
        << "-----------------------------------------" << nl
        << endl;

    // Main refinement engine
    meshRefinement meshRefiner
    (
        mesh,
        mergeDist,          // tolerance used in sorting coordinates
        overwrite,          // overwrite mesh files?
        surfaces,           // for surface intersection refinement
        features,           // for feature edges/point based refinement
        shells              // for volume (inside/outside) refinement
    );
    Info<< "Calculated surface intersections in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;

    // Some stats
    meshRefiner.printMeshInfo(debugLevel, "Initial mesh");

    meshRefiner.write
    (
        meshRefinement::debugType(debugLevel&meshRefinement::OBJINTERSECTIONS),
        meshRefinement::writeType(0),
        mesh.time().path()/meshRefiner.timeName()
    );


    // Add all the surface regions as patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //- Global surface region to patch (non faceZone surface) or patches
    //  (faceZone surfaces)
    labelList globalToMasterPatch;
    labelList globalToSlavePatch;
    {
        Info<< nl
            << "Adding patches for surface regions" << nl
            << "----------------------------------" << nl
            << endl;

        // From global region number to mesh patch.
        globalToMasterPatch.setSize(surfaces.nRegions(), -1);
        globalToSlavePatch.setSize(surfaces.nRegions(), -1);

        Info<< setf(ios_base::left)
            << setw(6) << "Patch"
            << setw(20) << "Type"
            << setw(30) << "Region" << nl
            << setw(6) << "-----"
            << setw(20) << "----"
            << setw(30) << "------" << endl;

        const labelList& surfaceGeometry = surfaces.surfaces();
        const PtrList<dictionary>& surfacePatchInfo = surfaces.patchInfo();

        forAll(surfaceGeometry, surfI)
        {
            label geomI = surfaceGeometry[surfI];

            const wordList& regNames = allGeometry.regionNames()[geomI];

            Info<< surfaces.names()[surfI] << ':' << nl << nl;

            if (surfaces.surfZones()[surfI].faceZoneName().empty())
            {
                // 'Normal' surface
                forAll(regNames, i)
                {
                    label globalRegionI = surfaces.globalRegion(surfI, i);

                    label patchi;

                    if (surfacePatchInfo.set(globalRegionI))
                    {
                        patchi = meshRefiner.addMeshedPatch
                        (
                            regNames[i],
                            surfacePatchInfo[globalRegionI]
                        );
                    }
                    else
                    {
                        dictionary patchInfo;
                        patchInfo.set("type", wallPolyPatch::typeName);

                        patchi = meshRefiner.addMeshedPatch
                        (
                            regNames[i],
                            patchInfo
                        );
                    }

                    Info<< setf(ios_base::left)
                        << setw(6) << patchi
                        << setw(20) << mesh.boundaryMesh()[patchi].type()
                        << setw(30) << regNames[i] << nl;

                    globalToMasterPatch[globalRegionI] = patchi;
                    globalToSlavePatch[globalRegionI] = patchi;
                }
            }
            else
            {
                // Zoned surface
                forAll(regNames, i)
                {
                    label globalRegionI = surfaces.globalRegion(surfI, i);

                    // Add master side patch
                    {
                        label patchi;

                        if (surfacePatchInfo.set(globalRegionI))
                        {
                            patchi = meshRefiner.addMeshedPatch
                            (
                                regNames[i],
                                surfacePatchInfo[globalRegionI]
                            );
                        }
                        else
                        {
                            dictionary patchInfo;
                            patchInfo.set("type", wallPolyPatch::typeName);

                            patchi = meshRefiner.addMeshedPatch
                            (
                                regNames[i],
                                patchInfo
                            );
                        }

                        Info<< setf(ios_base::left)
                            << setw(6) << patchi
                            << setw(20) << mesh.boundaryMesh()[patchi].type()
                            << setw(30) << regNames[i] << nl;

                        globalToMasterPatch[globalRegionI] = patchi;
                    }
                    // Add slave side patch
                    {
                        const word slaveName = regNames[i] + "_slave";
                        label patchi;

                        if (surfacePatchInfo.set(globalRegionI))
                        {
                            patchi = meshRefiner.addMeshedPatch
                            (
                                slaveName,
                                surfacePatchInfo[globalRegionI]
                            );
                        }
                        else
                        {
                            dictionary patchInfo;
                            patchInfo.set("type", wallPolyPatch::typeName);

                            patchi = meshRefiner.addMeshedPatch
                            (
                                slaveName,
                                patchInfo
                            );
                        }

                        Info<< setf(ios_base::left)
                            << setw(6) << patchi
                            << setw(20) << mesh.boundaryMesh()[patchi].type()
                            << setw(30) << slaveName << nl;

                        globalToSlavePatch[globalRegionI] = patchi;
                    }
                }
            }

            Info<< nl;
        }
        Info<< "Added patches in = "
            << mesh.time().cpuTimeIncrement() << " s" << nl << endl;
    }


    // Parallel
    // ~~~~~~~~

    // Decomposition
    autoPtr<decompositionMethod> decomposerPtr
    (
        decompositionMethod::New
        (
            decomposeDict
        )
    );
    decompositionMethod& decomposer = decomposerPtr();

    if (Pstream::parRun() && !decomposer.parallelAware())
    {
        FatalErrorInFunction
            << "You have selected decomposition method "
            << decomposer.typeName
            << " which is not parallel aware." << endl
            << "Please select one that is (hierarchical, ptscotch)"
            << exit(FatalError);
    }

    // Mesh distribution engine (uses tolerance to reconstruct meshes)
    fvMeshDistribute distributor(mesh);





    // Now do the real work -refinement -snapping -layers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const Switch wantRefine(meshDict.lookup("castellatedMesh"));
    const Switch wantSnap(meshDict.lookup("snap"));
    const Switch wantLayers(meshDict.lookup("addLayers"));

    // Refinement parameters
    const refinementParameters refineParams(refineDict);

    // Snap parameters
    const snapParameters snapParams(snapDict);

    // Layer addition parameters
    const layerParameters layerParams(layerDict, mesh.boundaryMesh());


    if (wantRefine)
    {
        cpuTime timer;

        snappyRefineDriver refineDriver
        (
            meshRefiner,
            decomposer,
            distributor,
            globalToMasterPatch,
            globalToSlavePatch
        );


        if (!overwrite && !debugLevel)
        {
            const_cast<Time&>(mesh.time())++;
        }


        refineDriver.doRefine
        (
            refineDict,
            refineParams,
            snapParams,
            refineParams.handleSnapProblems(),
            motionDict
        );


        if (!keepPatches && !wantSnap && !wantLayers)
        {
            removeZeroSizedPatches(mesh);
        }

        writeMesh
        (
            "Refined mesh",
            meshRefiner,
            debugLevel,
            meshRefinement::writeLevel()
        );

        Info<< "Mesh refined in = "
            << timer.cpuTimeIncrement() << " s." << endl;
    }

    if (wantSnap)
    {
        cpuTime timer;

        snappySnapDriver snapDriver
        (
            meshRefiner,
            globalToMasterPatch,
            globalToSlavePatch
        );

        if (!overwrite && !debugLevel)
        {
            const_cast<Time&>(mesh.time())++;
        }

        // Use the resolveFeatureAngle from the refinement parameters
        scalar curvature = refineParams.curvature();
        scalar planarAngle = refineParams.planarAngle();

        snapDriver.doSnap
        (
            snapDict,
            motionDict,
            curvature,
            planarAngle,
            snapParams
        );

        if (!keepPatches && !wantLayers)
        {
            removeZeroSizedPatches(mesh);
        }

        writeMesh
        (
            "Snapped mesh",
            meshRefiner,
            debugLevel,
            meshRefinement::writeLevel()
        );

        Info<< "Mesh snapped in = "
            << timer.cpuTimeIncrement() << " s." << endl;
    }

    if (wantLayers)
    {
        cpuTime timer;

        snappyLayerDriver layerDriver
        (
            meshRefiner,
            globalToMasterPatch,
            globalToSlavePatch
        );

        // Use the maxLocalCells from the refinement parameters
        bool preBalance = returnReduce
        (
            (mesh.nCells() >= refineParams.maxLocalCells()),
            orOp<bool>()
        );


        if (!overwrite && !debugLevel)
        {
            const_cast<Time&>(mesh.time())++;
        }

        layerDriver.doLayers
        (
            layerDict,
            motionDict,
            layerParams,
            preBalance,
            decomposer,
            distributor
        );

        if (!keepPatches)
        {
            removeZeroSizedPatches(mesh);
        }

        writeMesh
        (
            "Layer mesh",
            meshRefiner,
            debugLevel,
            meshRefinement::writeLevel()
        );

        Info<< "Layers added in = "
            << timer.cpuTimeIncrement() << " s." << endl;
    }


    {
        // Check final mesh
        Info<< "Checking final mesh ..." << endl;
        faceSet wrongFaces(mesh, "wrongFaces", mesh.nFaces()/100);
        motionSmoother::checkMesh(false, mesh, motionDict, wrongFaces);
        const label nErrors = returnReduce
        (
            wrongFaces.size(),
            sumOp<label>()
        );

        if (nErrors > 0)
        {
            Info<< "Finished meshing with " << nErrors << " illegal faces"
                << " (concave, zero area or negative cell pyramid volume)"
                << endl;
            wrongFaces.write();
        }
        else
        {
            Info<< "Finished meshing without any errors" << endl;
        }
    }


    if (surfaceSimplify)
    {
        const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

        labelHashSet includePatches(bMesh.size());

        if (args.optionFound("patches"))
        {
            includePatches = bMesh.patchSet
            (
                wordReList(args.optionLookup("patches")())
            );
        }
        else
        {
            forAll(bMesh, patchi)
            {
                const polyPatch& patch = bMesh[patchi];

                if (!isA<processorPolyPatch>(patch))
                {
                    includePatches.insert(patchi);
                }
            }
        }

        fileName outFileName
        (
            args.optionLookupOrDefault<fileName>
            (
                "outFile",
                "constant/triSurface/simplifiedSurface.stl"
            )
        );

        extractSurface
        (
            mesh,
            runTime,
            includePatches,
            outFileName
        );

        pointIOField cellCentres
        (
            IOobject
            (
                "internalCellCentres",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh.cellCentres()
        );

        cellCentres.write();
    }


    Info<< "Finished meshing in = "
        << runTime.elapsedCpuTime() << " s." << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
