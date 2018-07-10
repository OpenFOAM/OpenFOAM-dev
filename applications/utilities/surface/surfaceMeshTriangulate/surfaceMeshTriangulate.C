/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    surfaceMeshTriangulate

Description
    Extracts surface from a polyMesh. Depending on output surface format
    triangulates faces.

    Region numbers on faces cannot be guaranteed to be the same as the patch
    indices.

    Optionally only triangulates named patches.

    If run in parallel the processor patches get filtered out by default and
    the mesh gets merged (based on topology).

\*---------------------------------------------------------------------------*/

#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "ListListOps.H"
#include "uindirectPrimitivePatch.H"
#include "globalMeshData.H"
#include "globalIndex.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "extract surface from a polyMesh"
    );
    argList::validArgs.append("output surface file");
    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "excludeProcPatches",
        "exclude processor patches"
    );
    argList::addOption
    (
        "patches",
        "(patch0 .. patchN)",
        "only triangulate selected patches (wildcards supported)"
    );
    argList::addOption
    (
        "faceZones",
        "(fz0 .. fzN)",
        "triangulate selected faceZones (wildcards supported)"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const fileName outFileName(args[1]);

    Info<< "Extracting surface from boundaryMesh ..."
        << endl << endl;

    const bool includeProcPatches =
       !(
            args.optionFound("excludeProcPatches")
         || Pstream::parRun()
        );

    if (includeProcPatches)
    {
        Info<< "Including all processor patches." << nl << endl;
    }
    else if (Pstream::parRun())
    {
        Info<< "Excluding all processor patches." << nl << endl;
    }

    Info<< "Reading mesh from time " << runTime.value() << endl;

    #include "createNamedPolyMesh.H"


    // Create local surface from:
    // - explicitly named patches only (-patches (at your option)
    // - all patches (default in sequential mode)
    // - all non-processor patches (default in parallel mode)
    // - all non-processor patches (sequential mode, -excludeProcPatches
    //   (at your option)

    // Construct table of patches to include.
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

            if (includeProcPatches || !isA<processorPolyPatch>(patch))
            {
                includePatches.insert(patchi);
            }
        }
    }



    const faceZoneMesh& fzm = mesh.faceZones();
    labelHashSet includeFaceZones(fzm.size());

    if (args.optionFound("faceZones"))
    {
        wordReList zoneNames(args.optionLookup("faceZones")());
        const wordList allZoneNames(fzm.names());
        forAll(zoneNames, i)
        {
            const wordRe& zoneName = zoneNames[i];

            labelList zoneIDs = findStrings(zoneName, allZoneNames);

            forAll(zoneIDs, j)
            {
                includeFaceZones.insert(zoneIDs[j]);
            }

            if (zoneIDs.empty())
            {
                WarningInFunction
                    << "Cannot find any faceZone name matching "
                    << zoneName << endl;
            }

        }
        Info<< "Additionally triangulating faceZones "
            << UIndirectList<word>(allZoneNames, includeFaceZones.sortedToc())
            << endl;
    }



    // From (name of) patch to compact 'zone' index
    HashTable<label> compactZoneID(1000);
    // Mesh face and compact zone indx
    DynamicList<label> faceLabels;
    DynamicList<label> compactZones;

    {
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

        HashTable<label> zoneSize(1000);
        forAllConstIter(labelHashSet, includeFaceZones, iter)
        {
            const faceZone& pp = fzm[iter.key()];
            zoneSize.insert(pp.name(), pp.size());
            nFaces += pp.size();
        }


        Pstream::mapCombineGather(patchSize, plusEqOp<label>());
        Pstream::mapCombineGather(zoneSize, plusEqOp<label>());


        // Allocate compact numbering for all patches/faceZones
        forAllConstIter(HashTable<label>, patchSize, iter)
        {
            label sz = compactZoneID.size();
            compactZoneID.insert(iter.key(), sz);
        }

        forAllConstIter(HashTable<label>, zoneSize, iter)
        {
            label sz = compactZoneID.size();
            // Info<< "For faceZone " << iter.key() << " allocating zoneID "
            //    << sz << endl;
            compactZoneID.insert(iter.key(), sz);
        }


        Pstream::mapCombineScatter(compactZoneID);


        // Rework HashTable into labelList just for speed of conversion
        labelList patchToCompactZone(bMesh.size(), -1);
        labelList faceZoneToCompactZone(bMesh.size(), -1);
        forAllConstIter(HashTable<label>, compactZoneID, iter)
        {
            label patchi = bMesh.findPatchID(iter.key());
            if (patchi != -1)
            {
                patchToCompactZone[patchi] = iter();
            }
            else
            {
                label zoneI = fzm.findZoneID(iter.key());
                faceZoneToCompactZone[zoneI] = iter();
            }
        }


        faceLabels.setCapacity(nFaces);
        compactZones.setCapacity(nFaces);

        // Collect faces on patches
        forAllConstIter(labelHashSet, includePatches, iter)
        {
            const polyPatch& pp = bMesh[iter.key()];
            forAll(pp, i)
            {
                faceLabels.append(pp.start()+i);
                compactZones.append(patchToCompactZone[pp.index()]);
            }
        }
        // Collect faces on faceZones
        forAllConstIter(labelHashSet, includeFaceZones, iter)
        {
            const faceZone& pp = fzm[iter.key()];
            forAll(pp, i)
            {
                faceLabels.append(pp[i]);
                compactZones.append(faceZoneToCompactZone[pp.index()]);
            }
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
    gatheredZones[Pstream::myProcNo()] = compactZones.xfer();
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
            xferMove(allPoints),
            xferMove(allFaces),
            xferMove(allZones),
            xferMove(surfZones)
        );


        MeshedSurface<face> sortedFace(unsortedFace);

        fileName globalCasePath
        (
            runTime.processorCase()
          ? runTime.path()/".."/outFileName
          : runTime.path()/outFileName
        );
        globalCasePath.clean();

        Info<< "Writing merged surface to " << globalCasePath << endl;

        sortedFace.write(globalCasePath);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
