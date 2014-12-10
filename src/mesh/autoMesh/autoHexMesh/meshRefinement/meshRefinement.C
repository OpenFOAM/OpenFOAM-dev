/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "meshRefinement.H"
#include "volMesh.H"
#include "volFields.H"
#include "surfaceMesh.H"
#include "syncTools.H"
#include "Time.H"
#include "refinementHistory.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "decompositionMethod.H"
#include "regionSplit.H"
#include "fvMeshDistribute.H"
#include "indirectPrimitivePatch.H"
#include "polyTopoChange.H"
#include "removeCells.H"
#include "mapDistributePolyMesh.H"
#include "localPointRegion.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "slipPointPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "calculatedPointPatchFields.H"
#include "cyclicSlipPointPatchFields.H"
#include "processorPointPatch.H"
#include "globalIndex.H"
#include "meshTools.H"
#include "OFstream.H"
#include "geomDecomp.H"
#include "Random.H"
#include "searchableSurfaces.H"
#include "treeBoundBox.H"
#include "zeroGradientFvPatchFields.H"
#include "fvMeshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshRefinement, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::meshRefinement::IOdebugType,
        5
    >::names[] =
    {
        "mesh",
        //"scalarLevels",
        "intersections",
        "featureSeeds",
        "attraction",
        "layerInfo"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::meshRefinement::IOoutputType,
        1
    >::names[] =
    {
        "layerInfo"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::meshRefinement::IOwriteType,
        4
    >::names[] =
    {
        "mesh",
        "scalarLevels",
        "layerSets",
        "layerFields"
    };

}

const Foam::NamedEnum<Foam::meshRefinement::IOdebugType, 5>
Foam::meshRefinement::IOdebugTypeNames;

const Foam::NamedEnum<Foam::meshRefinement::IOoutputType, 1>
Foam::meshRefinement::IOoutputTypeNames;

const Foam::NamedEnum<Foam::meshRefinement::IOwriteType, 4>
Foam::meshRefinement::IOwriteTypeNames;


Foam::meshRefinement::writeType Foam::meshRefinement::writeLevel_;

Foam::meshRefinement::outputType Foam::meshRefinement::outputLevel_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshRefinement::calcNeighbourData
(
    labelList& neiLevel,
    pointField& neiCc
)  const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    label nBoundaryFaces = mesh_.nFaces() - mesh_.nInternalFaces();

    if (neiLevel.size() != nBoundaryFaces || neiCc.size() != nBoundaryFaces)
    {
        FatalErrorIn("meshRefinement::calcNeighbour(..)") << "nBoundaries:"
            << nBoundaryFaces << " neiLevel:" << neiLevel.size()
            << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    labelHashSet addedPatchIDSet(meshedPatches());

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        const labelUList& faceCells = pp.faceCells();
        const vectorField::subField faceCentres = pp.faceCentres();
        const vectorField::subField faceAreas = pp.faceAreas();

        label bFaceI = pp.start()-mesh_.nInternalFaces();

        if (pp.coupled())
        {
            forAll(faceCells, i)
            {
                neiLevel[bFaceI] = cellLevel[faceCells[i]];
                neiCc[bFaceI] = cellCentres[faceCells[i]];
                bFaceI++;
            }
        }
        else if (addedPatchIDSet.found(patchI))
        {
            // Face was introduced from cell-cell intersection. Try to
            // reconstruct other side cell(centre). Three possibilities:
            // - cells same size.
            // - preserved cell smaller. Not handled.
            // - preserved cell larger.
            forAll(faceCells, i)
            {
                // Extrapolate the face centre.
                vector fn = faceAreas[i];
                fn /= mag(fn)+VSMALL;

                label own = faceCells[i];
                label ownLevel = cellLevel[own];
                label faceLevel = meshCutter_.getAnchorLevel(pp.start()+i);

                // Normal distance from face centre to cell centre
                scalar d = ((faceCentres[i] - cellCentres[own]) & fn);
                if (faceLevel > ownLevel)
                {
                    // Other cell more refined. Adjust normal distance
                    d *= 0.5;
                }
                neiLevel[bFaceI] = faceLevel;
                // Calculate other cell centre by extrapolation
                neiCc[bFaceI] = faceCentres[i] + d*fn;
                bFaceI++;
            }
        }
        else
        {
            forAll(faceCells, i)
            {
                neiLevel[bFaceI] = cellLevel[faceCells[i]];
                neiCc[bFaceI] = faceCentres[i];
                bFaceI++;
            }
        }
    }

    // Swap coupled boundaries. Apply separation to cc since is coordinate.
    syncTools::swapBoundaryFacePositions(mesh_, neiCc);
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);
}


// Find intersections of edges (given by their two endpoints) with surfaces.
// Returns first intersection if there are more than one.
void Foam::meshRefinement::updateIntersections(const labelList& changedFaces)
{
    const pointField& cellCentres = mesh_.cellCentres();

    // Stats on edges to test. Count proc faces only once.
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));

    {
        label nMasterFaces = 0;
        forAll(isMasterFace, faceI)
        {
            if (isMasterFace.get(faceI) == 1)
            {
                nMasterFaces++;
            }
        }
        reduce(nMasterFaces, sumOp<label>());

        label nChangedFaces = 0;
        forAll(changedFaces, i)
        {
            if (isMasterFace.get(changedFaces[i]) == 1)
            {
                nChangedFaces++;
            }
        }
        reduce(nChangedFaces, sumOp<label>());

        Info<< "Edge intersection testing:" << nl
            << "    Number of edges             : " << nMasterFaces << nl
            << "    Number of edges to retest   : " << nChangedFaces
            << endl;
    }


    // Get boundary face centre and level. Coupled aware.
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(neiLevel, neiCc);

    // Collect segments we want to test for
    pointField start(changedFaces.size());
    pointField end(changedFaces.size());

    forAll(changedFaces, i)
    {
        label faceI = changedFaces[i];
        label own = mesh_.faceOwner()[faceI];

        start[i] = cellCentres[own];
        if (mesh_.isInternalFace(faceI))
        {
            end[i] = cellCentres[mesh_.faceNeighbour()[faceI]];
        }
        else
        {
            end[i] = neiCc[faceI-mesh_.nInternalFaces()];
        }
    }

    // Extend segments a bit
    {
        const vectorField smallVec(Foam::sqrt(SMALL)*(end-start));
        start -= smallVec;
        end += smallVec;
    }


    // Do tests in one go
    labelList surfaceHit;
    {
        labelList surfaceLevel;
        surfaces_.findHigherIntersection
        (
            start,
            end,
            labelList(start.size(), -1),    // accept any intersection
            surfaceHit,
            surfaceLevel
        );
    }

    // Keep just surface hit
    forAll(surfaceHit, i)
    {
        surfaceIndex_[changedFaces[i]] = surfaceHit[i];
    }

    // Make sure both sides have same information. This should be
    // case in general since same vectors but just to make sure.
    syncTools::syncFaceList(mesh_, surfaceIndex_, maxEqOp<label>());

    label nHits = countHits();
    label nTotHits = returnReduce(nHits, sumOp<label>());

    Info<< "    Number of intersected edges : " << nTotHits << endl;

    // Set files to same time as mesh
    setInstance(mesh_.facesInstance());
}


void Foam::meshRefinement::testSyncPointList
(
    const string& msg,
    const polyMesh& mesh,
    const List<scalar>& fld
)
{
    if (fld.size() != mesh.nPoints())
    {
        FatalErrorIn
        (
            "meshRefinement::testSyncPointList(const polyMesh&"
            ", const List<scalar>&)"
        )   << msg << endl
            << "fld size:" << fld.size() << " mesh points:" << mesh.nPoints()
            << abort(FatalError);
    }

    Pout<< "Checking field " << msg << endl;
    scalarField minFld(fld);
    syncTools::syncPointList
    (
        mesh,
        minFld,
        minEqOp<scalar>(),
        GREAT
    );
    scalarField maxFld(fld);
    syncTools::syncPointList
    (
        mesh,
        maxFld,
        maxEqOp<scalar>(),
        -GREAT
    );
    forAll(minFld, pointI)
    {
        const scalar& minVal = minFld[pointI];
        const scalar& maxVal = maxFld[pointI];
        if (mag(minVal-maxVal) > SMALL)
        {
            Pout<< msg << " at:" << mesh.points()[pointI] << nl
                << "    minFld:" << minVal << nl
                << "    maxFld:" << maxVal << nl
                << endl;
        }
    }
}


void Foam::meshRefinement::testSyncPointList
(
    const string& msg,
    const polyMesh& mesh,
    const List<point>& fld
)
{
    if (fld.size() != mesh.nPoints())
    {
        FatalErrorIn
        (
            "meshRefinement::testSyncPointList(const polyMesh&"
            ", const List<point>&)"
        )   << msg << endl
            << "fld size:" << fld.size() << " mesh points:" << mesh.nPoints()
            << abort(FatalError);
    }

    Pout<< "Checking field " << msg << endl;
    pointField minFld(fld);
    syncTools::syncPointList
    (
        mesh,
        minFld,
        minMagSqrEqOp<point>(),
        point(GREAT, GREAT, GREAT)
    );
    pointField maxFld(fld);
    syncTools::syncPointList
    (
        mesh,
        maxFld,
        maxMagSqrEqOp<point>(),
        vector::zero
    );
    forAll(minFld, pointI)
    {
        const point& minVal = minFld[pointI];
        const point& maxVal = maxFld[pointI];
        if (mag(minVal-maxVal) > SMALL)
        {
            Pout<< msg << " at:" << mesh.points()[pointI] << nl
                << "    minFld:" << minVal << nl
                << "    maxFld:" << maxVal << nl
                << endl;
        }
    }
}


void Foam::meshRefinement::checkData()
{
    Pout<< "meshRefinement::checkData() : Checking refinement structure."
        << endl;
    meshCutter_.checkMesh();

    Pout<< "meshRefinement::checkData() : Checking refinement levels."
        << endl;
    meshCutter_.checkRefinementLevels(1, labelList(0));


    label nBnd = mesh_.nFaces()-mesh_.nInternalFaces();

    Pout<< "meshRefinement::checkData() : Checking synchronization."
        << endl;

    // Check face centres
    {
        // Boundary face centres
        pointField::subList boundaryFc
        (
            mesh_.faceCentres(),
            nBnd,
            mesh_.nInternalFaces()
        );

        // Get neighbouring face centres
        pointField neiBoundaryFc(boundaryFc);
        syncTools::syncBoundaryFacePositions
        (
            mesh_,
            neiBoundaryFc,
            eqOp<point>()
        );

        // Compare
        testSyncBoundaryFaceList
        (
            mergeDistance_,
            "testing faceCentres : ",
            boundaryFc,
            neiBoundaryFc
        );
    }
    // Check meshRefinement
    {
        // Get boundary face centre and level. Coupled aware.
        labelList neiLevel(nBnd);
        pointField neiCc(nBnd);
        calcNeighbourData(neiLevel, neiCc);

        // Collect segments we want to test for
        pointField start(mesh_.nFaces());
        pointField end(mesh_.nFaces());

        forAll(start, faceI)
        {
            start[faceI] = mesh_.cellCentres()[mesh_.faceOwner()[faceI]];

            if (mesh_.isInternalFace(faceI))
            {
                end[faceI] = mesh_.cellCentres()[mesh_.faceNeighbour()[faceI]];
            }
            else
            {
                end[faceI] = neiCc[faceI-mesh_.nInternalFaces()];
            }
        }

        // Extend segments a bit
        {
            const vectorField smallVec(Foam::sqrt(SMALL)*(end-start));
            start -= smallVec;
            end += smallVec;
        }


        // Do tests in one go
        labelList surfaceHit;
        {
            labelList surfaceLevel;
            surfaces_.findHigherIntersection
            (
                start,
                end,
                labelList(start.size(), -1),    // accept any intersection
                surfaceHit,
                surfaceLevel
            );
        }
        // Get the coupled hit
        labelList neiHit
        (
            SubList<label>
            (
                surfaceHit,
                nBnd,
                mesh_.nInternalFaces()
            )
        );
        syncTools::swapBoundaryFaceList(mesh_, neiHit);

        // Check
        forAll(surfaceHit, faceI)
        {
            if (surfaceIndex_[faceI] != surfaceHit[faceI])
            {
                if (mesh_.isInternalFace(faceI))
                {
                    WarningIn("meshRefinement::checkData()")
                        << "Internal face:" << faceI
                        << " fc:" << mesh_.faceCentres()[faceI]
                        << " cached surfaceIndex_:" << surfaceIndex_[faceI]
                        << " current:" << surfaceHit[faceI]
                        << " ownCc:"
                        << mesh_.cellCentres()[mesh_.faceOwner()[faceI]]
                        << " neiCc:"
                        << mesh_.cellCentres()[mesh_.faceNeighbour()[faceI]]
                        << endl;
                }
                else if
                (
                    surfaceIndex_[faceI]
                 != neiHit[faceI-mesh_.nInternalFaces()]
                )
                {
                    WarningIn("meshRefinement::checkData()")
                        << "Boundary face:" << faceI
                        << " fc:" << mesh_.faceCentres()[faceI]
                        << " cached surfaceIndex_:" << surfaceIndex_[faceI]
                        << " current:" << surfaceHit[faceI]
                        << " ownCc:"
                        << mesh_.cellCentres()[mesh_.faceOwner()[faceI]]
                        << " end:" << end[faceI]
                        << endl;
                }
            }
        }
    }
    {
        labelList::subList boundarySurface
        (
            surfaceIndex_,
            mesh_.nFaces()-mesh_.nInternalFaces(),
            mesh_.nInternalFaces()
        );

        labelList neiBoundarySurface(boundarySurface);
        syncTools::swapBoundaryFaceList
        (
            mesh_,
            neiBoundarySurface
        );

        // Compare
        testSyncBoundaryFaceList
        (
            0,                              // tolerance
            "testing surfaceIndex() : ",
            boundarySurface,
            neiBoundarySurface
        );
    }


    // Find duplicate faces
    Pout<< "meshRefinement::checkData() : Counting duplicate faces."
        << endl;

    labelList duplicateFace
    (
        localPointRegion::findDuplicateFaces
        (
            mesh_,
            identity(mesh_.nFaces()-mesh_.nInternalFaces())
          + mesh_.nInternalFaces()
        )
    );

    // Count
    {
        label nDup = 0;

        forAll(duplicateFace, i)
        {
            if (duplicateFace[i] != -1)
            {
                nDup++;
            }
        }
        nDup /= 2;  // will have counted both faces of duplicate
        Pout<< "meshRefinement::checkData() : Found " << nDup
            << " duplicate pairs of faces." << endl;
    }
}


void Foam::meshRefinement::setInstance(const fileName& inst)
{
    meshCutter_.setInstance(inst);
    surfaceIndex_.instance() = inst;
}


// Remove cells. Put exposedFaces (output of getExposedFaces(cellsToRemove))
// into exposedPatchIDs.
Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::doRemoveCells
(
    const labelList& cellsToRemove,
    const labelList& exposedFaces,
    const labelList& exposedPatchIDs,
    removeCells& cellRemover
)
{
    polyTopoChange meshMod(mesh_);

    // Arbitrary: put exposed faces into last patch.
    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        exposedPatchIDs,
        meshMod
    );

    // Change the mesh (no inflation)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh (since morphing might not do this)
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes. No other way to do this?
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());
    setInstance(mesh_.facesInstance());

    // Update local mesh data
    cellRemover.updateMesh(map);

    // Update intersections. Recalculate intersections for exposed faces.
    labelList newExposedFaces = renumber
    (
        map().reverseFaceMap(),
        exposedFaces
    );

    //Pout<< "removeCells : updating intersections for "
    //    << newExposedFaces.size() << " newly exposed faces." << endl;

    updateMesh(map, newExposedFaces);

    return map;
}


// Split faces
Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::splitFaces
(
    const labelList& splitFaces,
    const labelPairList& splits
)
{
    polyTopoChange meshMod(mesh_);

    forAll(splitFaces, i)
    {
        label faceI = splitFaces[i];
        const face& f = mesh_.faces()[faceI];


        // Split as start and end index in face
        const labelPair& split = splits[i];

        label nVerts = split[1]-split[0];
        if (nVerts < 0)
        {
            nVerts += f.size();
        }
        nVerts += 1;


        // Split into f0, f1
        face f0(nVerts);

        label fp = split[0];
        forAll(f0, i)
        {
            f0[i] = f[fp];
            fp = f.fcIndex(fp);
        }

        face f1(f.size()-f0.size()+2);
        fp = split[1];
        forAll(f1, i)
        {
            f1[i] = f[fp];
            fp = f.fcIndex(fp);
        }


        // Determine face properties
        label own = mesh_.faceOwner()[faceI];
        label nei = -1;
        label patchI = -1;
        if (faceI >= mesh_.nInternalFaces())
        {
            patchI = mesh_.boundaryMesh().whichPatch(faceI);
        }
        else
        {
            nei = mesh_.faceNeighbour()[faceI];
        }

        label zoneI = mesh_.faceZones().whichZone(faceI);
        bool zoneFlip = false;
        if (zoneI != -1)
        {
            const faceZone& fz = mesh_.faceZones()[zoneI];
            zoneFlip = fz.flipMap()[fz.whichFace(faceI)];
        }


Pout<< "face:" << faceI << " verts:" << f
    << " split into f0:" << f0
    << " f1:" << f1 << endl;

        // Change/add faces
        meshMod.modifyFace
        (
            f0,                         // modified face
            faceI,                      // label of face
            own,                        // owner
            nei,                        // neighbour
            false,                      // face flip
            patchI,                     // patch for face
            zoneI,                      // zone for face
            zoneFlip                    // face flip in zone
        );
        meshMod.addFace
        (
            f1,                         // modified face
            own,                        // owner
            nei,                        // neighbour
            -1,                         // master point
            -1,                         // master edge
            faceI,                      // master face
            false,                      // face flip
            patchI,                     // patch for face
            zoneI,                      // zone for face
            zoneFlip                    // face flip in zone
        );
    }



    // Change the mesh (no inflation)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh (since morphing might not do this)
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes. No other way to do this?
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());
    setInstance(mesh_.facesInstance());

    // Update local mesh data
    const labelList& oldToNew = map().reverseFaceMap();
    labelList newSplitFaces(renumber(oldToNew, splitFaces));
    // Add added faces (every splitFaces becomes two faces)
    label sz = newSplitFaces.size();
    newSplitFaces.setSize(2*sz);
    forAll(map().faceMap(), faceI)
    {
        label oldFaceI = map().faceMap()[faceI];
        if (oldToNew[oldFaceI] != faceI)
        {
            // Added face
            newSplitFaces[sz++] = faceI;
        }
    }

    updateMesh(map, newSplitFaces);

    return map;
}


//// Determine for multi-processor regions the lowest numbered cell on the
//// lowest numbered processor.
//void Foam::meshRefinement::getCoupledRegionMaster
//(
//    const globalIndex& globalCells,
//    const boolList& blockedFace,
//    const regionSplit& globalRegion,
//    Map<label>& regionToMaster
//) const
//{
//    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
//
//    forAll(patches, patchI)
//    {
//        const polyPatch& pp = patches[patchI];
//
//        if (isA<processorPolyPatch>(pp))
//        {
//            forAll(pp, i)
//            {
//                label faceI = pp.start()+i;
//
//                if (!blockedFace[faceI])
//                {
//                    // Only if there is a connection to the neighbour
//                    // will there be a multi-domain region; if not through
//                    // this face then through another.
//
//                    label cellI = mesh_.faceOwner()[faceI];
//                    label globalCellI = globalCells.toGlobal(cellI);
//
//                    Map<label>::iterator iter =
//                        regionToMaster.find(globalRegion[cellI]);
//
//                    if (iter != regionToMaster.end())
//                    {
//                        label master = iter();
//                        iter() = min(master, globalCellI);
//                    }
//                    else
//                    {
//                        regionToMaster.insert
//                        (
//                            globalRegion[cellI],
//                            globalCellI
//                        );
//                    }
//                }
//            }
//        }
//    }
//
//    // Do reduction
//    Pstream::mapCombineGather(regionToMaster, minEqOp<label>());
//    Pstream::mapCombineScatter(regionToMaster);
//}
//
//
//void Foam::meshRefinement::calcLocalRegions
//(
//    const globalIndex& globalCells,
//    const labelList& globalRegion,
//    const Map<label>& coupledRegionToMaster,
//    const scalarField& cellWeights,
//
//    Map<label>& globalToLocalRegion,
//    pointField& localPoints,
//    scalarField& localWeights
//) const
//{
//    globalToLocalRegion.resize(globalRegion.size());
//    DynamicList<point> localCc(globalRegion.size()/2);
//    DynamicList<scalar> localWts(globalRegion.size()/2);
//
//    forAll(globalRegion, cellI)
//    {
//        Map<label>::const_iterator fndMaster =
//            coupledRegionToMaster.find(globalRegion[cellI]);
//
//        if (fndMaster != coupledRegionToMaster.end())
//        {
//            // Multi-processor region.
//            if (globalCells.toGlobal(cellI) == fndMaster())
//            {
//                // I am master. Allocate region for me.
//                globalToLocalRegion.insert
//                (
//                    globalRegion[cellI],
//                    localCc.size()
//                );
//                localCc.append(mesh_.cellCentres()[cellI]);
//                localWts.append(cellWeights[cellI]);
//            }
//        }
//        else
//        {
//            // Single processor region.
//            if
//            (
//                globalToLocalRegion.insert
//                (
//                    globalRegion[cellI],
//                    localCc.size()
//                )
//            )
//            {
//                localCc.append(mesh_.cellCentres()[cellI]);
//                localWts.append(cellWeights[cellI]);
//            }
//        }
//    }
//
//    localPoints.transfer(localCc);
//    localWeights.transfer(localWts);
//
//    if (localPoints.size() != globalToLocalRegion.size())
//    {
//        FatalErrorIn("calcLocalRegions(..)")
//            << "localPoints:" << localPoints.size()
//            << " globalToLocalRegion:" << globalToLocalRegion.size()
//            << exit(FatalError);
//    }
//}
//
//
//Foam::label Foam::meshRefinement::getShiftedRegion
//(
//    const globalIndex& indexer,
//    const Map<label>& globalToLocalRegion,
//    const Map<label>& coupledRegionToShifted,
//    const label globalRegion
//)
//{
//    Map<label>::const_iterator iter =
//        globalToLocalRegion.find(globalRegion);
//
//    if (iter != globalToLocalRegion.end())
//    {
//        // Region is 'owned' locally. Convert local region index into global.
//        return indexer.toGlobal(iter());
//    }
//    else
//    {
//        return coupledRegionToShifted[globalRegion];
//    }
//}
//
//
//// Add if not yet present
//void Foam::meshRefinement::addUnique(const label elem, labelList& lst)
//{
//    if (findIndex(lst, elem) == -1)
//    {
//        label sz = lst.size();
//        lst.setSize(sz+1);
//        lst[sz] = elem;
//    }
//}
//
//
//void Foam::meshRefinement::calcRegionRegions
//(
//    const labelList& globalRegion,
//    const Map<label>& globalToLocalRegion,
//    const Map<label>& coupledRegionToMaster,
//    labelListList& regionRegions
//) const
//{
//    // Global region indexing since we now know the shifted regions.
//    globalIndex shiftIndexer(globalToLocalRegion.size());
//
//    // Redo the coupledRegionToMaster to be in shifted region indexing.
//    Map<label> coupledRegionToShifted(coupledRegionToMaster.size());
//    forAllConstIter(Map<label>, coupledRegionToMaster, iter)
//    {
//        label region = iter.key();
//
//        Map<label>::const_iterator fndRegion = globalToLocalRegion.find
//        (region);
//
//        if (fndRegion != globalToLocalRegion.end())
//        {
//            // A local cell is master of this region. Get its shifted region.
//            coupledRegionToShifted.insert
//            (
//                region,
//                shiftIndexer.toGlobal(fndRegion())
//            );
//        }
//        Pstream::mapCombineGather(coupledRegionToShifted, minEqOp<label>());
//        Pstream::mapCombineScatter(coupledRegionToShifted);
//    }
//
//
//    // Determine region-region connectivity.
//    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    // This is for all my regions (so my local ones or the ones I am master
//    // of) the neighbouring region indices.
//
//
//    // Transfer lists.
//    PtrList<HashSet<edge, Hash<edge> > > regionConnectivity
//    (Pstream::nProcs());
//    forAll(regionConnectivity, procI)
//    {
//        if (procI != Pstream::myProcNo())
//        {
//            regionConnectivity.set
//            (
//                procI,
//                new HashSet<edge, Hash<edge> >
//                (
//                    coupledRegionToShifted.size()
//                  / Pstream::nProcs()
//                )
//            );
//        }
//    }
//
//
//    // Connectivity. For all my local regions the connected regions.
//    regionRegions.setSize(globalToLocalRegion.size());
//
//    // Add all local connectivity to regionRegions; add all non-local
//    // to the transferlists.
//    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
//    {
//        label ownRegion = globalRegion[mesh_.faceOwner()[faceI]];
//        label neiRegion = globalRegion[mesh_.faceNeighbour()[faceI]];
//
//        if (ownRegion != neiRegion)
//        {
//            label shiftOwnRegion = getShiftedRegion
//            (
//                shiftIndexer,
//                globalToLocalRegion,
//                coupledRegionToShifted,
//                ownRegion
//            );
//            label shiftNeiRegion = getShiftedRegion
//            (
//                shiftIndexer,
//                globalToLocalRegion,
//                coupledRegionToShifted,
//                neiRegion
//            );
//
//
//            // Connection between two regions. Send to owner of region.
//            // - is ownRegion 'owned' by me
//            // - is neiRegion 'owned' by me
//
//            if (shiftIndexer.isLocal(shiftOwnRegion))
//            {
//                label localI = shiftIndexer.toLocal(shiftOwnRegion);
//                addUnique(shiftNeiRegion, regionRegions[localI]);
//            }
//            else
//            {
//                label masterProc = shiftIndexer.whichProcID(shiftOwnRegion);
//                regionConnectivity[masterProc].insert
//                (
//                    edge(shiftOwnRegion, shiftNeiRegion)
//                );
//            }
//
//            if (shiftIndexer.isLocal(shiftNeiRegion))
//            {
//                label localI = shiftIndexer.toLocal(shiftNeiRegion);
//                addUnique(shiftOwnRegion, regionRegions[localI]);
//            }
//            else
//            {
//                label masterProc = shiftIndexer.whichProcID(shiftNeiRegion);
//                regionConnectivity[masterProc].insert
//                (
//                    edge(shiftOwnRegion, shiftNeiRegion)
//                );
//            }
//        }
//    }
//
//
//    // Send
//    forAll(regionConnectivity, procI)
//    {
//        if (procI != Pstream::myProcNo())
//        {
//            OPstream str(Pstream::blocking, procI);
//            str << regionConnectivity[procI];
//        }
//    }
//    // Receive
//    forAll(regionConnectivity, procI)
//    {
//        if (procI != Pstream::myProcNo())
//        {
//            IPstream str(Pstream::blocking, procI);
//            str >> regionConnectivity[procI];
//        }
//    }
//
//    // Add to addressing.
//    forAll(regionConnectivity, procI)
//    {
//        if (procI != Pstream::myProcNo())
//        {
//            for
//            (
//                HashSet<edge, Hash<edge> >::const_iterator iter =
//                    regionConnectivity[procI].begin();
//                iter != regionConnectivity[procI].end();
//                ++iter
//            )
//            {
//                const edge& e = iter.key();
//
//                bool someLocal = false;
//                if (shiftIndexer.isLocal(e[0]))
//                {
//                    label localI = shiftIndexer.toLocal(e[0]);
//                    addUnique(e[1], regionRegions[localI]);
//                    someLocal = true;
//                }
//                if (shiftIndexer.isLocal(e[1]))
//                {
//                    label localI = shiftIndexer.toLocal(e[1]);
//                    addUnique(e[0], regionRegions[localI]);
//                    someLocal = true;
//                }
//
//                if (!someLocal)
//                {
//                    FatalErrorIn("calcRegionRegions(..)")
//                        << "Received from processor " << procI
//                        << " connection " << e
//                        << " where none of the elements is local to me."
//                        << abort(FatalError);
//                }
//            }
//        }
//    }
//}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::meshRefinement::meshRefinement
(
    fvMesh& mesh,
    const scalar mergeDistance,
    const bool overwrite,
    const refinementSurfaces& surfaces,
    const refinementFeatures& features,
    const shellSurfaces& shells
)
:
    mesh_(mesh),
    mergeDistance_(mergeDistance),
    overwrite_(overwrite),
    oldInstance_(mesh.pointsInstance()),
    surfaces_(surfaces),
    features_(features),
    shells_(shells),
    meshCutter_
    (
        mesh,
        false   // do not try to read history.
    ),
    surfaceIndex_
    (
        IOobject
        (
            "surfaceIndex",
            mesh_.facesInstance(),
            fvMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        labelList(mesh_.nFaces(), -1)
    ),
    userFaceData_(0)
{
    // recalculate intersections for all faces
    updateIntersections(identity(mesh_.nFaces()));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::meshRefinement::countHits() const
{
    // Stats on edges to test. Count proc faces only once.
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));

    label nHits = 0;

    forAll(surfaceIndex_, faceI)
    {
        if (surfaceIndex_[faceI] >= 0 && isMasterFace.get(faceI) == 1)
        {
            nHits++;
        }
    }
    return nHits;
}


//// Determine distribution to move connected regions onto one processor.
//Foam::labelList Foam::meshRefinement::decomposeCombineRegions
//(
//    const scalarField& cellWeights,
//    const boolList& blockedFace,
//    const List<labelPair>& explicitConnections,
//    decompositionMethod& decomposer
//) const
//{
//    // Determine global regions, separated by blockedFaces
//    regionSplit globalRegion(mesh_, blockedFace, explicitConnections);
//
//    // Now globalRegion has global region per cell. Problem is that
//    // the region might span multiple domains so we want to get
//    // a "region master" per domain. Note that multi-processor
//    // regions can only occur on cells on coupled patches.
//    // Note: since the number of regions does not change by this the
//    // process can be seen as just a shift of a region onto a single
//    // processor.
//
//
//    // Global cell numbering engine
//    globalIndex globalCells(mesh_.nCells());
//
//
//    // Determine per coupled region the master cell (lowest numbered cell
//    // on lowest numbered processor)
//    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    // (does not determine master for non-coupled=fully-local regions)
//
//    Map<label> coupledRegionToMaster(mesh_.nFaces()-mesh_.nInternalFaces());
//    getCoupledRegionMaster
//    (
//        globalCells,
//        blockedFace,
//        globalRegion,
//        coupledRegionToMaster
//    );
//
//    // Determine my regions
//    // ~~~~~~~~~~~~~~~~~~~~
//    // These are the fully local ones or the coupled ones of which I am master
//
//    Map<label> globalToLocalRegion;
//    pointField localPoints;
//    scalarField localWeights;
//    calcLocalRegions
//    (
//        globalCells,
//        globalRegion,
//        coupledRegionToMaster,
//        cellWeights,
//
//        globalToLocalRegion,
//        localPoints,
//        localWeights
//    );
//
//
//
//    // Find distribution for regions
//    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//    labelList regionDistribution;
//
//    if (isA<geomDecomp>(decomposer))
//    {
//        regionDistribution = decomposer.decompose(localPoints, localWeights);
//    }
//    else
//    {
//        labelListList regionRegions;
//        calcRegionRegions
//        (
//            globalRegion,
//            globalToLocalRegion,
//            coupledRegionToMaster,
//
//            regionRegions
//        );
//
//        regionDistribution = decomposer.decompose
//        (
//            regionRegions,
//            localPoints,
//            localWeights
//        );
//    }
//
//
//
//    // Convert region-based decomposition back to cell-based one
//    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//    // Transfer destination processor back to all. Use global reduce for now.
//    Map<label> regionToDist(coupledRegionToMaster.size());
//    forAllConstIter(Map<label>, coupledRegionToMaster, iter)
//    {
//        label region = iter.key();
//
//        Map<label>::const_iterator regionFnd = globalToLocalRegion.find
//        (region);
//
//        if (regionFnd != globalToLocalRegion.end())
//        {
//            // Master cell is local. Store my distribution.
//            regionToDist.insert(iter.key(), regionDistribution[regionFnd()]);
//        }
//        else
//        {
//            // Master cell is not on this processor. Make sure it is
//            // overridden.
//            regionToDist.insert(iter.key(), labelMax);
//        }
//    }
//    Pstream::mapCombineGather(regionToDist, minEqOp<label>());
//    Pstream::mapCombineScatter(regionToDist);
//
//
//    // Determine destination for all cells
//    labelList distribution(mesh_.nCells());
//
//    forAll(globalRegion, cellI)
//    {
//        Map<label>::const_iterator fndRegion =
//            regionToDist.find(globalRegion[cellI]);
//
//        if (fndRegion != regionToDist.end())
//        {
//            distribution[cellI] = fndRegion();
//        }
//        else
//        {
//            // region is local to the processor.
//            label localRegionI = globalToLocalRegion[globalRegion[cellI]];
//
//            distribution[cellI] = regionDistribution[localRegionI];
//        }
//    }
//
//    return distribution;
//}


Foam::autoPtr<Foam::mapDistributePolyMesh> Foam::meshRefinement::balance
(
    const bool keepZoneFaces,
    const bool keepBaffles,
    const scalarField& cellWeights,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor
)
{
    autoPtr<mapDistributePolyMesh> map;

    if (Pstream::parRun())
    {
        //if (debug_)
        //{
        //    const_cast<Time&>(mesh_.time())++;
        //}

        // Wanted distribution
        labelList distribution;


        // Faces where owner and neighbour are not 'connected' so can
        // go to different processors.
        boolList blockedFace;
        label nUnblocked = 0;

        // Faces that move as block onto single processor
        PtrList<labelList> specifiedProcessorFaces;
        labelList specifiedProcessor;

        // Pairs of baffles
        List<labelPair> couples;

        // Constraints from decomposeParDict
        decomposer.setConstraints
        (
            mesh_,
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            couples
        );


        if (keepZoneFaces || keepBaffles)
        {
            if (keepZoneFaces)
            {
                // Determine decomposition to keep/move surface zones
                // on one processor. The reason is that snapping will make these
                // into baffles, move and convert them back so if they were
                // proc boundaries after baffling&moving the points might be no
                // longer synchronised so recoupling will fail. To prevent this
                // keep owner&neighbour of such a surface zone on the same
                // processor.

                const PtrList<surfaceZonesInfo>& surfZones =
                    surfaces().surfZones();
                const faceZoneMesh& fZones = mesh_.faceZones();
                const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

                // Get faces whose owner and neighbour should stay together,
                // i.e. they are not 'blocked'.

                forAll(surfZones, surfI)
                {
                    const word& fzName = surfZones[surfI].faceZoneName();

                    if (fzName.size())
                    {
                        // Get zone
                        const faceZone& fZone = fZones[fzName];

                        forAll(fZone, i)
                        {
                            label faceI = fZone[i];
                            if (blockedFace[faceI])
                            {
                                if
                                (
                                    mesh_.isInternalFace(faceI)
                                 || pbm[pbm.whichPatch(faceI)].coupled()
                                )
                                {
                                    blockedFace[faceI] = false;
                                    nUnblocked++;
                                }
                            }
                        }
                    }
                }


                // If the faceZones are not synchronised the blockedFace
                // might not be synchronised. If you are sure the faceZones
                // are synchronised remove below check.
                syncTools::syncFaceList
                (
                    mesh_,
                    blockedFace,
                    andEqOp<bool>()     // combine operator
                );
            }
            reduce(nUnblocked, sumOp<label>());

            if (keepZoneFaces)
            {
                Info<< "Found " << nUnblocked
                    << " zoned faces to keep together." << endl;
            }


            // Extend un-blockedFaces with any cyclics
            {
                boolList separatedCoupledFace(mesh_.nFaces(), false);
                selectSeparatedCoupledFaces(separatedCoupledFace);

                label nSeparated = 0;
                forAll(separatedCoupledFace, faceI)
                {
                    if (separatedCoupledFace[faceI])
                    {
                        if (blockedFace[faceI])
                        {
                            blockedFace[faceI] = false;
                            nSeparated++;
                        }
                    }
                }
                reduce(nSeparated, sumOp<label>());
                Info<< "Found " << nSeparated
                    << " separated coupled faces to keep together." << endl;

                nUnblocked += nSeparated;
            }


            if (keepBaffles)
            {
                label nBnd = mesh_.nFaces()-mesh_.nInternalFaces();

                labelList coupledFace(mesh_.nFaces(), -1);
                {
                    // Get boundary baffles that need to stay together
                    List<labelPair> allCouples
                    (
                        localPointRegion::findDuplicateFacePairs(mesh_)
                    );

                    // Merge with any couples from
                    // decompositionMethod::setConstraints
                    forAll(couples, i)
                    {
                        const labelPair& baffle = couples[i];
                        coupledFace[baffle.first()] = baffle.second();
                        coupledFace[baffle.second()] = baffle.first();
                    }
                    forAll(allCouples, i)
                    {
                        const labelPair& baffle = allCouples[i];
                        coupledFace[baffle.first()] = baffle.second();
                        coupledFace[baffle.second()] = baffle.first();
                    }
                }

                couples.setSize(nBnd);
                label nCpl = 0;
                forAll(coupledFace, faceI)
                {
                    if (coupledFace[faceI] != -1 && faceI < coupledFace[faceI])
                    {
                        couples[nCpl++] = labelPair(faceI, coupledFace[faceI]);
                    }
                }
                couples.setSize(nCpl);
            }
            label nCouples = returnReduce(couples.size(), sumOp<label>());

            if (keepBaffles)
            {
                Info<< "Found " << nCouples << " baffles to keep together."
                    << endl;
            }

            //if (nUnblocked > 0 || nCouples > 0)
            //{
            //    Info<< "Applying special decomposition to keep baffles"
            //        << " and zoned faces together." << endl;
            //
            //    distribution = decomposeCombineRegions
            //    (
            //        cellWeights,
            //        blockedFace,
            //        couples,
            //        decomposer
            //    );
            //
            //    labelList nProcCells(distributor.countCells(distribution));
            //    Pstream::listCombineGather(nProcCells, plusEqOp<label>());
            //    Pstream::listCombineScatter(nProcCells);
            //
            //    Info<< "Calculated decomposition:" << endl;
            //    forAll(nProcCells, procI)
            //    {
            //        Info<< "    " << procI << '\t' << nProcCells[procI]
            //            << endl;
            //    }
            //    Info<< endl;
            //}
            //else
            //{
            //    // Normal decomposition
            //    distribution = decomposer.decompose
            //    (
            //        mesh_,
            //        mesh_.cellCentres(),
            //        cellWeights
            //    );
            //}
        }
        //else
        //{
        //    // Normal decomposition
        //    distribution = decomposer.decompose
        //    (
        //        mesh_,
        //        cellWeights
        //    );
        //}


        // Make sure blockedFace not set on couples
        forAll(couples, i)
        {
            const labelPair& baffle = couples[i];
            blockedFace[baffle.first()] = false;
            blockedFace[baffle.second()] = false;
        }

        distribution = decomposer.decompose
        (
            mesh_,
            cellWeights,
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            couples                 // explicit connections
        );

        if (debug)
        {
            labelList nProcCells(distributor.countCells(distribution));
            Pout<< "Wanted distribution:" << nProcCells << endl;

            Pstream::listCombineGather(nProcCells, plusEqOp<label>());
            Pstream::listCombineScatter(nProcCells);

            Pout<< "Wanted resulting decomposition:" << endl;
            forAll(nProcCells, procI)
            {
                Pout<< "    " << procI << '\t' << nProcCells[procI] << endl;
            }
            Pout<< endl;
        }
        // Do actual sending/receiving of mesh
        map = distributor.distribute(distribution);

        // Update numbering of meshRefiner
        distribute(map);

        // Set correct instance (for if overwrite)
        mesh_.setInstance(timeName());
        setInstance(mesh_.facesInstance());
    }

    return map;
}


// Helper function to get intersected faces
Foam::labelList Foam::meshRefinement::intersectedFaces() const
{
    label nBoundaryFaces = 0;

    forAll(surfaceIndex_, faceI)
    {
        if (surfaceIndex_[faceI] != -1)
        {
            nBoundaryFaces++;
        }
    }

    labelList surfaceFaces(nBoundaryFaces);
    nBoundaryFaces = 0;

    forAll(surfaceIndex_, faceI)
    {
        if (surfaceIndex_[faceI] != -1)
        {
            surfaceFaces[nBoundaryFaces++] = faceI;
        }
    }
    return surfaceFaces;
}


// Helper function to get points used by faces
Foam::labelList Foam::meshRefinement::intersectedPoints() const
{
    const faceList& faces = mesh_.faces();

    // Mark all points on faces that will become baffles
    PackedBoolList isBoundaryPoint(mesh_.nPoints());
    label nBoundaryPoints = 0;

    forAll(surfaceIndex_, faceI)
    {
        if (surfaceIndex_[faceI] != -1)
        {
            const face& f = faces[faceI];

            forAll(f, fp)
            {
                if (isBoundaryPoint.set(f[fp], 1u))
                {
                    nBoundaryPoints++;
                }
            }
        }
    }

    //// Insert all meshed patches.
    //labelList adaptPatchIDs(meshedPatches());
    //forAll(adaptPatchIDs, i)
    //{
    //    label patchI = adaptPatchIDs[i];
    //
    //    if (patchI != -1)
    //    {
    //        const polyPatch& pp = mesh_.boundaryMesh()[patchI];
    //
    //        label faceI = pp.start();
    //
    //        forAll(pp, i)
    //        {
    //            const face& f = faces[faceI];
    //
    //            forAll(f, fp)
    //            {
    //                if (isBoundaryPoint.set(f[fp], 1u))
    //                    nBoundaryPoints++;
    //                }
    //            }
    //            faceI++;
    //        }
    //    }
    //}


    // Pack
    labelList boundaryPoints(nBoundaryPoints);
    nBoundaryPoints = 0;
    forAll(isBoundaryPoint, pointI)
    {
        if (isBoundaryPoint.get(pointI) == 1u)
        {
            boundaryPoints[nBoundaryPoints++] = pointI;
        }
    }

    return boundaryPoints;
}


//- Create patch from set of patches
Foam::autoPtr<Foam::indirectPrimitivePatch> Foam::meshRefinement::makePatch
(
    const polyMesh& mesh,
    const labelList& patchIDs
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Count faces.
    label nFaces = 0;

    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];

        nFaces += pp.size();
    }

    // Collect faces.
    labelList addressing(nFaces);
    nFaces = 0;

    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];

        label meshFaceI = pp.start();

        forAll(pp, i)
        {
            addressing[nFaces++] = meshFaceI++;
        }
    }

    return autoPtr<indirectPrimitivePatch>
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>(mesh.faces(), addressing),
            mesh.points()
        )
    );
}


// Construct pointVectorField with correct boundary conditions
Foam::tmp<Foam::pointVectorField> Foam::meshRefinement::makeDisplacementField
(
    const pointMesh& pMesh,
    const labelList& adaptPatchIDs
)
{
    const polyMesh& mesh = pMesh();

    // Construct displacement field.
    const pointBoundaryMesh& pointPatches = pMesh.boundary();

    wordList patchFieldTypes
    (
        pointPatches.size(),
        slipPointPatchVectorField::typeName
    );

    forAll(adaptPatchIDs, i)
    {
        patchFieldTypes[adaptPatchIDs[i]] =
            fixedValuePointPatchVectorField::typeName;
    }

    forAll(pointPatches, patchI)
    {
        if (isA<processorPointPatch>(pointPatches[patchI]))
        {
            patchFieldTypes[patchI] = calculatedPointPatchVectorField::typeName;
        }
        else if (isA<cyclicPointPatch>(pointPatches[patchI]))
        {
            patchFieldTypes[patchI] = cyclicSlipPointPatchVectorField::typeName;
        }
    }

    // Note: time().timeName() instead of meshRefinement::timeName() since
    // postprocessable field.
    tmp<pointVectorField> tfld
    (
        new pointVectorField
        (
            IOobject
            (
                "pointDisplacement",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedVector("displacement", dimLength, vector::zero),
            patchFieldTypes
        )
    );
    return tfld;
}


void Foam::meshRefinement::checkCoupledFaceZones(const polyMesh& mesh)
{
    const faceZoneMesh& fZones = mesh.faceZones();

    // Check any zones are present anywhere and in same order

    {
        List<wordList> zoneNames(Pstream::nProcs());
        zoneNames[Pstream::myProcNo()] = fZones.names();
        Pstream::gatherList(zoneNames);
        Pstream::scatterList(zoneNames);
        // All have same data now. Check.
        forAll(zoneNames, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                if (zoneNames[procI] != zoneNames[Pstream::myProcNo()])
                {
                    FatalErrorIn
                    (
                        "meshRefinement::checkCoupledFaceZones(const polyMesh&)"
                    )   << "faceZones are not synchronised on processors." << nl
                        << "Processor " << procI << " has faceZones "
                        << zoneNames[procI] << nl
                        << "Processor " << Pstream::myProcNo()
                        << " has faceZones "
                        << zoneNames[Pstream::myProcNo()] << nl
                        << exit(FatalError);
                }
            }
        }
    }

    // Check that coupled faces are present on both sides.

    labelList faceToZone(mesh.nFaces()-mesh.nInternalFaces(), -1);

    forAll(fZones, zoneI)
    {
        const faceZone& fZone = fZones[zoneI];

        forAll(fZone, i)
        {
            label bFaceI = fZone[i]-mesh.nInternalFaces();

            if (bFaceI >= 0)
            {
                if (faceToZone[bFaceI] == -1)
                {
                    faceToZone[bFaceI] = zoneI;
                }
                else if (faceToZone[bFaceI] == zoneI)
                {
                    FatalErrorIn
                    (
                        "meshRefinement::checkCoupledFaceZones(const polyMesh&)"
                    )   << "Face " << fZone[i] << " in zone "
                        << fZone.name()
                        << " is twice in zone!"
                        << abort(FatalError);
                }
                else
                {
                    FatalErrorIn
                    (
                        "meshRefinement::checkCoupledFaceZones(const polyMesh&)"
                    )   << "Face " << fZone[i] << " in zone "
                        << fZone.name()
                        << " is also in zone "
                        << fZones[faceToZone[bFaceI]].name()
                        << abort(FatalError);
                }
            }
        }
    }

    labelList neiFaceToZone(faceToZone);
    syncTools::swapBoundaryFaceList(mesh, neiFaceToZone);

    forAll(faceToZone, i)
    {
        if (faceToZone[i] != neiFaceToZone[i])
        {
            FatalErrorIn
            (
                "meshRefinement::checkCoupledFaceZones(const polyMesh&)"
            )   << "Face " << mesh.nInternalFaces()+i
                << " is in zone " << faceToZone[i]
                << ", its coupled face is in zone " << neiFaceToZone[i]
                << abort(FatalError);
        }
    }
}


void Foam::meshRefinement::calculateEdgeWeights
(
    const polyMesh& mesh,
    const PackedBoolList& isMasterEdge,
    const labelList& meshPoints,
    const edgeList& edges,
    scalarField& edgeWeights,
    scalarField& invSumWeight
)
{
    const pointField& pts = mesh.points();

    // Calculate edgeWeights and inverse sum of edge weights
    edgeWeights.setSize(isMasterEdge.size());
    invSumWeight.setSize(meshPoints.size());

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];
        scalar eMag = max
        (
            SMALL,
            mag
            (
                pts[meshPoints[e[1]]]
              - pts[meshPoints[e[0]]]
            )
        );
        edgeWeights[edgeI] = 1.0/eMag;
    }

    // Sum per point all edge weights
    weightedSum
    (
        mesh,
        isMasterEdge,
        meshPoints,
        edges,
        edgeWeights,
        scalarField(meshPoints.size(), 1.0),  // data
        invSumWeight
    );

    // Inplace invert
    forAll(invSumWeight, pointI)
    {
        scalar w = invSumWeight[pointI];

        if (w > 0.0)
        {
            invSumWeight[pointI] = 1.0/w;
        }
    }
}


Foam::label Foam::meshRefinement::appendPatch
(
    fvMesh& mesh,
    const label insertPatchI,
    const word& patchName,
    const dictionary& patchDict
)
{
    // Clear local fields and e.g. polyMesh parallelInfo.
    mesh.clearOut();

    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    label patchI = polyPatches.size();

    // Add polyPatch at the end
    polyPatches.setSize(patchI+1);
    polyPatches.set
    (
        patchI,
        polyPatch::New
        (
            patchName,
            patchDict,
            insertPatchI,
            polyPatches
        )
    );
    fvPatches.setSize(patchI+1);
    fvPatches.set
    (
        patchI,
        fvPatch::New
        (
            polyPatches[patchI],  // point to newly added polyPatch
            mesh.boundary()
        )
    );

    addPatchFields<volScalarField>
    (
        mesh,
        calculatedFvPatchField<scalar>::typeName
    );
    addPatchFields<volVectorField>
    (
        mesh,
        calculatedFvPatchField<vector>::typeName
    );
    addPatchFields<volSphericalTensorField>
    (
        mesh,
        calculatedFvPatchField<sphericalTensor>::typeName
    );
    addPatchFields<volSymmTensorField>
    (
        mesh,
        calculatedFvPatchField<symmTensor>::typeName
    );
    addPatchFields<volTensorField>
    (
        mesh,
        calculatedFvPatchField<tensor>::typeName
    );

    // Surface fields

    addPatchFields<surfaceScalarField>
    (
        mesh,
        calculatedFvPatchField<scalar>::typeName
    );
    addPatchFields<surfaceVectorField>
    (
        mesh,
        calculatedFvPatchField<vector>::typeName
    );
    addPatchFields<surfaceSphericalTensorField>
    (
        mesh,
        calculatedFvPatchField<sphericalTensor>::typeName
    );
    addPatchFields<surfaceSymmTensorField>
    (
        mesh,
        calculatedFvPatchField<symmTensor>::typeName
    );
    addPatchFields<surfaceTensorField>
    (
        mesh,
        calculatedFvPatchField<tensor>::typeName
    );
    return patchI;
}


// Adds patch if not yet there. Returns patchID.
Foam::label Foam::meshRefinement::addPatch
(
    fvMesh& mesh,
    const word& patchName,
    const dictionary& patchInfo
)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    const label patchI = polyPatches.findPatchID(patchName);
    if (patchI != -1)
    {
        // Already there
        return patchI;
    }


    label insertPatchI = polyPatches.size();
    label startFaceI = mesh.nFaces();

    forAll(polyPatches, patchI)
    {
        const polyPatch& pp = polyPatches[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            insertPatchI = patchI;
            startFaceI = pp.start();
            break;
        }
    }

    dictionary patchDict(patchInfo);
    patchDict.set("nFaces", 0);
    patchDict.set("startFace", startFaceI);

    // Below is all quite a hack. Feel free to change once there is a better
    // mechanism to insert and reorder patches.

    label addedPatchI = appendPatch(mesh, insertPatchI, patchName, patchDict);


    // Create reordering list
    // patches before insert position stay as is
    labelList oldToNew(addedPatchI+1);
    for (label i = 0; i < insertPatchI; i++)
    {
        oldToNew[i] = i;
    }
    // patches after insert position move one up
    for (label i = insertPatchI; i < addedPatchI; i++)
    {
        oldToNew[i] = i+1;
    }
    // appended patch gets moved to insert position
    oldToNew[addedPatchI] = insertPatchI;

    // Shuffle into place
    polyPatches.reorder(oldToNew, true);
    fvPatches.reorder(oldToNew);

    reorderPatchFields<volScalarField>(mesh, oldToNew);
    reorderPatchFields<volVectorField>(mesh, oldToNew);
    reorderPatchFields<volSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<volSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<volTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceScalarField>(mesh, oldToNew);
    reorderPatchFields<surfaceVectorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceTensorField>(mesh, oldToNew);

    return insertPatchI;
}


Foam::label Foam::meshRefinement::addMeshedPatch
(
    const word& name,
    const dictionary& patchInfo
)
{
    label meshedI = findIndex(meshedPatches_, name);

    if (meshedI != -1)
    {
        // Already there. Get corresponding polypatch
        return mesh_.boundaryMesh().findPatchID(name);
    }
    else
    {
        // Add patch
        label patchI = addPatch(mesh_, name, patchInfo);

//        dictionary patchDict(patchInfo);
//        patchDict.set("nFaces", 0);
//        patchDict.set("startFace", 0);
//        autoPtr<polyPatch> ppPtr
//        (
//            polyPatch::New
//            (
//                name,
//                patchDict,
//                0,
//                mesh_.boundaryMesh()
//            )
//        );
//        label patchI = fvMeshTools::addPatch
//        (
//            mesh_,
//            ppPtr(),
//            dictionary(),       // optional field values
//            calculatedFvPatchField<scalar>::typeName,
//            true
//        );

        // Store
        meshedPatches_.append(name);

        return patchI;
    }
}


Foam::labelList Foam::meshRefinement::meshedPatches() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    DynamicList<label> patchIDs(meshedPatches_.size());
    forAll(meshedPatches_, i)
    {
        label patchI = patches.findPatchID(meshedPatches_[i]);

        if (patchI == -1)
        {
            FatalErrorIn("meshRefinement::meshedPatches() const")
                << "Problem : did not find patch " << meshedPatches_[i]
                << endl << "Valid patches are " << patches.names()
                << abort(FatalError);
        }
        if (!polyPatch::constraintType(patches[patchI].type()))
        {
            patchIDs.append(patchI);
        }
    }

    return patchIDs;
}


void Foam::meshRefinement::selectSeparatedCoupledFaces(boolList& selected) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        // Check all coupled. Avoid using .coupled() so we also pick up AMI.
        if (isA<coupledPolyPatch>(patches[patchI]))
        {
            const coupledPolyPatch& cpp = refCast<const coupledPolyPatch>
            (
                patches[patchI]
            );

            if (cpp.separated() || !cpp.parallel())
            {
                forAll(cpp, i)
                {
                    selected[cpp.start()+i] = true;
                }
            }
        }
    }
}


Foam::label Foam::meshRefinement::findRegion
(
    const polyMesh& mesh,
    const labelList& cellToRegion,
    const vector& perturbVec,
    const point& p
)
{
    label regionI = -1;
    label cellI = mesh.findCell(p);
    if (cellI != -1)
    {
        regionI = cellToRegion[cellI];
    }
    reduce(regionI, maxOp<label>());

    if (regionI == -1)
    {
        // See if we can perturb a bit
        cellI = mesh.findCell(p+perturbVec);
        if (cellI != -1)
        {
            regionI = cellToRegion[cellI];
        }
        reduce(regionI, maxOp<label>());
    }
    return regionI;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::splitMeshRegions
(
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const point& keepPoint
)
{
    // Force calculation of face decomposition (used in findCell)
    (void)mesh_.tetBasePtIs();

    // Determine connected regions. regionSplit is the labelList with the
    // region per cell.

    boolList blockedFace(mesh_.nFaces(), false);
    selectSeparatedCoupledFaces(blockedFace);

    regionSplit cellRegion(mesh_, blockedFace);

    label regionI = findRegion
    (
        mesh_,
        cellRegion,
        mergeDistance_*vector(1,1,1), // note:1,1,1 should really be normalised
        keepPoint
    );

    if (regionI == -1)
    {
        FatalErrorIn
        (
            "meshRefinement::splitMeshRegions(const point&)"
        )   << "Point " << keepPoint
            << " is not inside the mesh." << nl
            << "Bounding box of the mesh:" << mesh_.bounds()
            << exit(FatalError);
    }

    // Subset
    // ~~~~~~

    // Get cells to remove
    DynamicList<label> cellsToRemove(mesh_.nCells());
    forAll(cellRegion, cellI)
    {
        if (cellRegion[cellI] != regionI)
        {
            cellsToRemove.append(cellI);
        }
    }
    cellsToRemove.shrink();

    label nCellsToKeep = mesh_.nCells() - cellsToRemove.size();
    reduce(nCellsToKeep, sumOp<label>());

    Info<< "Keeping all cells in region " << regionI
        << " containing point " << keepPoint << endl
        << "Selected for keeping : "
        << nCellsToKeep
        << " cells." << endl;


    // Remove cells
    removeCells cellRemover(mesh_);

    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    labelList exposedPatch;

    label nExposedFaces = returnReduce(exposedFaces.size(), sumOp<label>());
    if (nExposedFaces)
    {
        //FatalErrorIn
        //(
        //    "meshRefinement::splitMeshRegions(const point&)"
        //)   << "Removing non-reachable cells should only expose"
        //    << " boundary faces" << nl
        //    << "ExposedFaces:" << exposedFaces << abort(FatalError);

        // Patch for exposed faces for lack of anything sensible.
        label defaultPatch = 0;
        if (globalToMasterPatch.size())
        {
            defaultPatch = globalToMasterPatch[0];
        }

        WarningIn
        (
            "meshRefinement::splitMeshRegions(const point&)"
        )   << "Removing non-reachable cells exposes "
            << nExposedFaces << " internal or coupled faces." << endl
            << "    These get put into patch " << defaultPatch << endl;
        exposedPatch.setSize(exposedFaces.size(), defaultPatch);
    }

    return doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatch,
        cellRemover
    );
}


void Foam::meshRefinement::distribute(const mapDistributePolyMesh& map)
{
    // mesh_ already distributed; distribute my member data

    // surfaceQueries_ ok.

    // refinement
    meshCutter_.distribute(map);

    // surfaceIndex is face data.
    map.distributeFaceData(surfaceIndex_);

    // maintainedFaces are indices of faces.
    forAll(userFaceData_, i)
    {
        map.distributeFaceData(userFaceData_[i].second());
    }

    // Redistribute surface and any fields on it.
    {
        Random rndGen(653213);

        // Get local mesh bounding box. Single box for now.
        List<treeBoundBox> meshBb(1);
        treeBoundBox& bb = meshBb[0];
        bb = treeBoundBox(mesh_.points());
        bb = bb.extend(rndGen, 1e-4);

        // Distribute all geometry (so refinementSurfaces and shellSurfaces)
        searchableSurfaces& geometry =
            const_cast<searchableSurfaces&>(surfaces_.geometry());

        forAll(geometry, i)
        {
            autoPtr<mapDistribute> faceMap;
            autoPtr<mapDistribute> pointMap;
            geometry[i].distribute
            (
                meshBb,
                false,          // do not keep outside triangles
                faceMap,
                pointMap
            );

            if (faceMap.valid())
            {
                // (ab)use the instance() to signal current modification time
                geometry[i].instance() = geometry[i].time().timeName();
            }

            faceMap.clear();
            pointMap.clear();
        }
    }
}


// Update local data for a mesh change
void Foam::meshRefinement::updateMesh
(
    const mapPolyMesh& map,
    const labelList& changedFaces
)
{
    Map<label> dummyMap(0);

    updateMesh(map, changedFaces, dummyMap, dummyMap, dummyMap);
}


void Foam::meshRefinement::storeData
(
    const labelList& pointsToStore,
    const labelList& facesToStore,
    const labelList& cellsToStore
)
{
    // For now only meshCutter has storable/retrievable data.
    meshCutter_.storeData
    (
        pointsToStore,
        facesToStore,
        cellsToStore
    );
}


void Foam::meshRefinement::updateMesh
(
    const mapPolyMesh& map,
    const labelList& changedFaces,
    const Map<label>& pointsToRestore,
    const Map<label>& facesToRestore,
    const Map<label>& cellsToRestore
)
{
    // For now only meshCutter has storable/retrievable data.

    // Update numbering of cells/vertices.
    meshCutter_.updateMesh
    (
        map,
        pointsToRestore,
        facesToRestore,
        cellsToRestore
    );

    // Update surfaceIndex
    updateList(map.faceMap(), label(-1), surfaceIndex_);

    // Update cached intersection information
    updateIntersections(changedFaces);

    // Update maintained faces
    forAll(userFaceData_, i)
    {
        labelList& data = userFaceData_[i].second();

        if (userFaceData_[i].first() == KEEPALL)
        {
            // extend list with face-from-face data
            updateList(map.faceMap(), label(-1), data);
        }
        else if (userFaceData_[i].first() == MASTERONLY)
        {
            // keep master only
            labelList newFaceData(map.faceMap().size(), -1);

            forAll(newFaceData, faceI)
            {
                label oldFaceI = map.faceMap()[faceI];

                if (oldFaceI >= 0 && map.reverseFaceMap()[oldFaceI] == faceI)
                {
                    newFaceData[faceI] = data[oldFaceI];
                }
            }
            data.transfer(newFaceData);
        }
        else
        {
            // remove any face that has been refined i.e. referenced more than
            // once.

            // 1. Determine all old faces that get referenced more than once.
            // These get marked with -1 in reverseFaceMap
            labelList reverseFaceMap(map.reverseFaceMap());

            forAll(map.faceMap(), faceI)
            {
                label oldFaceI = map.faceMap()[faceI];

                if (oldFaceI >= 0)
                {
                    if (reverseFaceMap[oldFaceI] != faceI)
                    {
                        // faceI is slave face. Mark old face.
                        reverseFaceMap[oldFaceI] = -1;
                    }
                }
            }

            // 2. Map only faces with intact reverseFaceMap
            labelList newFaceData(map.faceMap().size(), -1);
            forAll(newFaceData, faceI)
            {
                label oldFaceI = map.faceMap()[faceI];

                if (oldFaceI >= 0)
                {
                    if (reverseFaceMap[oldFaceI] == faceI)
                    {
                        newFaceData[faceI] = data[oldFaceI];
                    }
                }
            }
            data.transfer(newFaceData);
        }
    }
}


bool Foam::meshRefinement::write() const
{
    bool writeOk =
        mesh_.write()
     && meshCutter_.write()
     && surfaceIndex_.write();


    // Make sure that any distributed surfaces (so ones which probably have
    // been changed) get written as well.
    // Note: should ideally have some 'modified' flag to say whether it
    // has been changed or not.
    searchableSurfaces& geometry =
        const_cast<searchableSurfaces&>(surfaces_.geometry());

    forAll(geometry, i)
    {
        searchableSurface& s = geometry[i];

        // Check if instance() of surface is not constant or system.
        // Is good hint that surface is distributed.
        if
        (
            s.instance() != s.time().system()
         && s.instance() != s.time().caseSystem()
         && s.instance() != s.time().constant()
         && s.instance() != s.time().caseConstant()
        )
        {
            // Make sure it gets written to current time, not constant.
            s.instance() = s.time().timeName();
            writeOk = writeOk && s.write();
        }
    }

    return writeOk;
}


Foam::PackedBoolList Foam::meshRefinement::getMasterPoints
(
    const polyMesh& mesh,
    const labelList& meshPoints
)
{
    const globalIndex globalPoints(meshPoints.size());

    labelList myPoints(meshPoints.size());
    forAll(meshPoints, pointI)
    {
        myPoints[pointI] = globalPoints.toGlobal(pointI);
    }

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        myPoints,
        minEqOp<label>(),
        labelMax
    );


    PackedBoolList isPatchMasterPoint(meshPoints.size());
    forAll(meshPoints, pointI)
    {
        if (myPoints[pointI] == globalPoints.toGlobal(pointI))
        {
            isPatchMasterPoint[pointI] = true;
        }
    }

    return isPatchMasterPoint;
}


Foam::PackedBoolList Foam::meshRefinement::getMasterEdges
(
    const polyMesh& mesh,
    const labelList& meshEdges
)
{
    const globalIndex globalEdges(meshEdges.size());

    labelList myEdges(meshEdges.size());
    forAll(meshEdges, edgeI)
    {
        myEdges[edgeI] = globalEdges.toGlobal(edgeI);
    }

    syncTools::syncEdgeList
    (
        mesh,
        meshEdges,
        myEdges,
        minEqOp<label>(),
        labelMax
    );


    PackedBoolList isMasterEdge(meshEdges.size());
    forAll(meshEdges, edgeI)
    {
        if (myEdges[edgeI] == globalEdges.toGlobal(edgeI))
        {
            isMasterEdge[edgeI] = true;
        }
    }

    return isMasterEdge;
}


void Foam::meshRefinement::printMeshInfo(const bool debug, const string& msg)
const
{
    const globalMeshData& pData = mesh_.globalData();

    if (debug)
    {
        Pout<< msg.c_str()
            << " : cells(local):" << mesh_.nCells()
            << "  faces(local):" << mesh_.nFaces()
            << "  points(local):" << mesh_.nPoints()
            << endl;
    }

    {
        PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));
        label nMasterFaces = 0;
        forAll(isMasterFace, i)
        {
            if (isMasterFace[i])
            {
                nMasterFaces++;
            }
        }

        PackedBoolList isMeshMasterPoint(syncTools::getMasterPoints(mesh_));
        label nMasterPoints = 0;
        forAll(isMeshMasterPoint, i)
        {
            if (isMeshMasterPoint[i])
            {
                nMasterPoints++;
            }
        }

        Info<< msg.c_str()
            << " : cells:" << pData.nTotalCells()
            << "  faces:" << returnReduce(nMasterFaces, sumOp<label>())
            << "  points:" << returnReduce(nMasterPoints, sumOp<label>())
            << endl;
    }


    //if (debug)
    {
        const labelList& cellLevel = meshCutter_.cellLevel();

        labelList nCells(gMax(cellLevel)+1, 0);

        forAll(cellLevel, cellI)
        {
            nCells[cellLevel[cellI]]++;
        }

        Pstream::listCombineGather(nCells, plusEqOp<label>());
        Pstream::listCombineScatter(nCells);

        Info<< "Cells per refinement level:" << endl;
        forAll(nCells, levelI)
        {
            Info<< "    " << levelI << '\t' << nCells[levelI]
                << endl;
        }
    }
}


//- Return either time().constant() or oldInstance
Foam::word Foam::meshRefinement::timeName() const
{
    if (overwrite_ && mesh_.time().timeIndex() == 0)
    {
        return oldInstance_;
    }
    else
    {
        return mesh_.time().timeName();
    }
}


void Foam::meshRefinement::dumpRefinementLevel() const
{
    // Note: use time().timeName(), not meshRefinement::timeName()
    // so as to dump the fields to 0, not to constant.
    {
        volScalarField volRefLevel
        (
            IOobject
            (
                "cellLevel",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );

        const labelList& cellLevel = meshCutter_.cellLevel();

        forAll(volRefLevel, cellI)
        {
            volRefLevel[cellI] = cellLevel[cellI];
        }

        volRefLevel.write();
    }

    // Dump pointLevel
    {
        const pointMesh& pMesh = pointMesh::New(mesh_);

        pointScalarField pointRefLevel
        (
            IOobject
            (
                "pointLevel",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pMesh,
            dimensionedScalar("zero", dimless, 0)
        );

        const labelList& pointLevel = meshCutter_.pointLevel();

        forAll(pointRefLevel, pointI)
        {
            pointRefLevel[pointI] = pointLevel[pointI];
        }

        pointRefLevel.write();
    }
}


void Foam::meshRefinement::dumpIntersections(const fileName& prefix) const
{
    {
        const pointField& cellCentres = mesh_.cellCentres();

        OFstream str(prefix + "_edges.obj");
        label vertI = 0;
        Pout<< "meshRefinement::dumpIntersections :"
            << " Writing cellcentre-cellcentre intersections to file "
            << str.name() << endl;


        // Redo all intersections
        // ~~~~~~~~~~~~~~~~~~~~~~

        // Get boundary face centre and level. Coupled aware.
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
        pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
        calcNeighbourData(neiLevel, neiCc);

        labelList intersectionFaces(intersectedFaces());

        // Collect segments we want to test for
        pointField start(intersectionFaces.size());
        pointField end(intersectionFaces.size());

        forAll(intersectionFaces, i)
        {
            label faceI = intersectionFaces[i];
            start[i] = cellCentres[mesh_.faceOwner()[faceI]];

            if (mesh_.isInternalFace(faceI))
            {
                end[i] = cellCentres[mesh_.faceNeighbour()[faceI]];
            }
            else
            {
                end[i] = neiCc[faceI-mesh_.nInternalFaces()];
            }
        }

        // Extend segments a bit
        {
            const vectorField smallVec(Foam::sqrt(SMALL)*(end-start));
            start -= smallVec;
            end += smallVec;
        }


        // Do tests in one go
        labelList surfaceHit;
        List<pointIndexHit> surfaceHitInfo;
        surfaces_.findAnyIntersection
        (
            start,
            end,
            surfaceHit,
            surfaceHitInfo
        );

        forAll(intersectionFaces, i)
        {
            if (surfaceHit[i] != -1)
            {
                meshTools::writeOBJ(str, start[i]);
                vertI++;
                meshTools::writeOBJ(str, surfaceHitInfo[i].hitPoint());
                vertI++;
                meshTools::writeOBJ(str, end[i]);
                vertI++;
                str << "l " << vertI-2 << ' ' << vertI-1 << nl
                    << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }

    Pout<< endl;
}


void Foam::meshRefinement::write
(
    const debugType debugFlags,
    const writeType writeFlags,
    const fileName& prefix
) const
{
    if (writeFlags & WRITEMESH)
    {
        write();
    }
    if (writeFlags & WRITELEVELS)
    {
        dumpRefinementLevel();
    }
    if (debugFlags & OBJINTERSECTIONS && prefix.size())
    {
        dumpIntersections(prefix);
    }
}


Foam::meshRefinement::writeType Foam::meshRefinement::writeLevel()
{
    return writeLevel_;
}


void Foam::meshRefinement::writeLevel(const writeType flags)
{
    writeLevel_ = flags;
}


Foam::meshRefinement::outputType Foam::meshRefinement::outputLevel()
{
    return outputLevel_;
}


void Foam::meshRefinement::outputLevel(const outputType flags)
{
    outputLevel_ = flags;
}


// ************************************************************************* //
