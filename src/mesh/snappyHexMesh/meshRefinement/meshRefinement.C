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
        5
    >::names[] =
    {
        "mesh",
        "noRefinement",
        "scalarLevels",
        "layerSets",
        "layerFields"
    };

}

const Foam::NamedEnum<Foam::meshRefinement::IOdebugType, 5>
Foam::meshRefinement::IOdebugTypeNames;

const Foam::NamedEnum<Foam::meshRefinement::IOoutputType, 1>
Foam::meshRefinement::IOoutputTypeNames;

const Foam::NamedEnum<Foam::meshRefinement::IOwriteType, 5>
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
        FatalErrorInFunction
            << nBoundaryFaces << " neiLevel:" << neiLevel.size()
            << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    labelHashSet addedPatchIDSet(meshedPatches());

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        const labelUList& faceCells = pp.faceCells();
        const vectorField::subField faceCentres = pp.faceCentres();
        const vectorField::subField faceAreas = pp.faceAreas();

        label bFacei = pp.start()-mesh_.nInternalFaces();

        if (pp.coupled())
        {
            forAll(faceCells, i)
            {
                neiLevel[bFacei] = cellLevel[faceCells[i]];
                neiCc[bFacei] = cellCentres[faceCells[i]];
                bFacei++;
            }
        }
        else if (addedPatchIDSet.found(patchi))
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
                fn /= mag(fn)+vSmall;

                label own = faceCells[i];
                label ownLevel = cellLevel[own];
                label faceLevel = meshCutter_.faceLevel(pp.start()+i);

                // Normal distance from face centre to cell centre
                scalar d = ((faceCentres[i] - cellCentres[own]) & fn);
                if (faceLevel > ownLevel)
                {
                    // Other cell more refined. Adjust normal distance
                    d *= 0.5;
                }
                neiLevel[bFacei] = faceLevel;
                // Calculate other cell centre by extrapolation
                neiCc[bFacei] = faceCentres[i] + d*fn;
                bFacei++;
            }
        }
        else
        {
            forAll(faceCells, i)
            {
                neiLevel[bFacei] = cellLevel[faceCells[i]];
                neiCc[bFacei] = faceCentres[i];
                bFacei++;
            }
        }
    }

    // Swap coupled boundaries. Apply separation to cc since is coordinate.
    syncTools::swapBoundaryFacePositions(mesh_, neiCc);
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);
}


void Foam::meshRefinement::updateIntersections(const labelList& changedFaces)
{
    const pointField& cellCentres = mesh_.cellCentres();

    // Stats on edges to test. Count proc faces only once.
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));

    {
        label nMasterFaces = 0;
        forAll(isMasterFace, facei)
        {
            if (isMasterFace.get(facei) == 1)
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
        label facei = changedFaces[i];
        label own = mesh_.faceOwner()[facei];

        start[i] = cellCentres[own];
        if (mesh_.isInternalFace(facei))
        {
            end[i] = cellCentres[mesh_.faceNeighbour()[facei]];
        }
        else
        {
            end[i] = neiCc[facei-mesh_.nInternalFaces()];
        }
    }

    // Extend segments a bit
    {
        const vectorField smallVec(rootSmall*(end-start));
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
        FatalErrorInFunction
            << msg << endl
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
        great
    );
    scalarField maxFld(fld);
    syncTools::syncPointList
    (
        mesh,
        maxFld,
        maxEqOp<scalar>(),
        -great
    );
    forAll(minFld, pointi)
    {
        const scalar& minVal = minFld[pointi];
        const scalar& maxVal = maxFld[pointi];
        if (mag(minVal-maxVal) > small)
        {
            Pout<< msg << " at:" << mesh.points()[pointi] << nl
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
        FatalErrorInFunction
            << msg << endl
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
        point(great, great, great)
    );
    pointField maxFld(fld);
    syncTools::syncPointList
    (
        mesh,
        maxFld,
        maxMagSqrEqOp<point>(),
        vector::zero
    );
    forAll(minFld, pointi)
    {
        const point& minVal = minFld[pointi];
        const point& maxVal = maxFld[pointi];
        if (mag(minVal-maxVal) > small)
        {
            Pout<< msg << " at:" << mesh.points()[pointi] << nl
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

        forAll(start, facei)
        {
            start[facei] = mesh_.cellCentres()[mesh_.faceOwner()[facei]];

            if (mesh_.isInternalFace(facei))
            {
                end[facei] = mesh_.cellCentres()[mesh_.faceNeighbour()[facei]];
            }
            else
            {
                end[facei] = neiCc[facei-mesh_.nInternalFaces()];
            }
        }

        // Extend segments a bit
        {
            const vectorField smallVec(rootSmall*(end-start));
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
        forAll(surfaceHit, facei)
        {
            if (surfaceIndex_[facei] != surfaceHit[facei])
            {
                if (mesh_.isInternalFace(facei))
                {
                    WarningInFunction
                        << "Internal face:" << facei
                        << " fc:" << mesh_.faceCentres()[facei]
                        << " cached surfaceIndex_:" << surfaceIndex_[facei]
                        << " current:" << surfaceHit[facei]
                        << " ownCc:"
                        << mesh_.cellCentres()[mesh_.faceOwner()[facei]]
                        << " neiCc:"
                        << mesh_.cellCentres()[mesh_.faceNeighbour()[facei]]
                        << endl;
                }
                else if
                (
                    surfaceIndex_[facei]
                 != neiHit[facei-mesh_.nInternalFaces()]
                )
                {
                    WarningInFunction
                        << "Boundary face:" << facei
                        << " fc:" << mesh_.faceCentres()[facei]
                        << " cached surfaceIndex_:" << surfaceIndex_[facei]
                        << " current:" << surfaceHit[facei]
                        << " ownCc:"
                        << mesh_.cellCentres()[mesh_.faceOwner()[facei]]
                        << " end:" << end[facei]
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

    // Pout<< "removeCells : updating intersections for "
    //    << newExposedFaces.size() << " newly exposed faces." << endl;

    updateMesh(map, newExposedFaces);

    return map;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::splitFaces
(
    const labelList& splitFaces,
    const labelPairList& splits
)
{
    polyTopoChange meshMod(mesh_);

    forAll(splitFaces, i)
    {
        label facei = splitFaces[i];
        const face& f = mesh_.faces()[facei];

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
        label own = mesh_.faceOwner()[facei];
        label nei = -1;
        label patchi = -1;
        if (facei >= mesh_.nInternalFaces())
        {
            patchi = mesh_.boundaryMesh().whichPatch(facei);
        }
        else
        {
            nei = mesh_.faceNeighbour()[facei];
        }

        label zoneI = mesh_.faceZones().whichZone(facei);
        bool zoneFlip = false;
        if (zoneI != -1)
        {
            const faceZone& fz = mesh_.faceZones()[zoneI];
            zoneFlip = fz.flipMap()[fz.whichFace(facei)];
        }


        if (debug)
        {
            Pout<< "face:" << facei << " verts:" << f
                << " split into f0:" << f0
                << " f1:" << f1 << endl;
        }

        // Change/add faces
        meshMod.modifyFace
        (
            f0,                         // modified face
            facei,                      // label of face
            own,                        // owner
            nei,                        // neighbour
            false,                      // face flip
            patchi,                     // patch for face
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
            facei,                      // master face
            false,                      // face flip
            patchi,                     // patch for face
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
    forAll(map().faceMap(), facei)
    {
        label oldFacei = map().faceMap()[facei];
        if (oldToNew[oldFacei] != facei)
        {
            // Added face
            newSplitFaces[sz++] = facei;
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
//    forAll(patches, patchi)
//    {
//        const polyPatch& pp = patches[patchi];
//
//        if (isA<processorPolyPatch>(pp))
//        {
//            forAll(pp, i)
//            {
//                label facei = pp.start()+i;
//
//                if (!blockedFace[facei])
//                {
//                    // Only if there is a connection to the neighbour
//                    // will there be a multi-domain region; if not through
//                    // this face then through another.
//
//                    label celli = mesh_.faceOwner()[facei];
//                    label globalCelli = globalCells.toGlobal(celli);
//
//                    Map<label>::iterator iter =
//                        regionToMaster.find(globalRegion[celli]);
//
//                    if (iter != regionToMaster.end())
//                    {
//                        label master = iter();
//                        iter() = min(master, globalCelli);
//                    }
//                    else
//                    {
//                        regionToMaster.insert
//                        (
//                            globalRegion[celli],
//                            globalCelli
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
//    forAll(globalRegion, celli)
//    {
//        Map<label>::const_iterator fndMaster =
//            coupledRegionToMaster.find(globalRegion[celli]);
//
//        if (fndMaster != coupledRegionToMaster.end())
//        {
//            // Multi-processor region.
//            if (globalCells.toGlobal(celli) == fndMaster())
//            {
//                // I am master. Allocate region for me.
//                globalToLocalRegion.insert
//                (
//                    globalRegion[celli],
//                    localCc.size()
//                );
//                localCc.append(mesh_.cellCentres()[celli]);
//                localWts.append(cellWeights[celli]);
//            }
//        }
//        else
//        {
//            // Single processor region.
//            if
//            (
//                globalToLocalRegion.insert
//                (
//                    globalRegion[celli],
//                    localCc.size()
//                )
//            )
//            {
//                localCc.append(mesh_.cellCentres()[celli]);
//                localWts.append(cellWeights[celli]);
//            }
//        }
//    }
//
//    localPoints.transfer(localCc);
//    localWeights.transfer(localWts);
//
//    if (localPoints.size() != globalToLocalRegion.size())
//    {
//        FatalErrorInFunction
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
//    PtrList<HashSet<edge, Hash<edge>>> regionConnectivity
//    (Pstream::nProcs());
//    forAll(regionConnectivity, proci)
//    {
//        if (proci != Pstream::myProcNo())
//        {
//            regionConnectivity.set
//            (
//                proci,
//                new HashSet<edge, Hash<edge>>
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
//    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
//    {
//        label ownRegion = globalRegion[mesh_.faceOwner()[facei]];
//        label neiRegion = globalRegion[mesh_.faceNeighbour()[facei]];
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
//    forAll(regionConnectivity, proci)
//    {
//        if (proci != Pstream::myProcNo())
//        {
//            OPstream str(Pstream::commsTypes::blocking, proci);
//            str << regionConnectivity[proci];
//        }
//    }
//    // Receive
//    forAll(regionConnectivity, proci)
//    {
//        if (proci != Pstream::myProcNo())
//        {
//            IPstream str(Pstream::commsTypes::blocking, proci);
//            str >> regionConnectivity[proci];
//        }
//    }
//
//    // Add to addressing.
//    forAll(regionConnectivity, proci)
//    {
//        if (proci != Pstream::myProcNo())
//        {
//            for
//            (
//                HashSet<edge, Hash<edge>>::const_iterator iter =
//                    regionConnectivity[proci].begin();
//                iter != regionConnectivity[proci].end();
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
//                    FatalErrorInFunction
//                        << "Received from processor " << proci
//                        << " connection " << e
//                        << " where none of the elements is local to me."
//                        << abort(FatalError);
//                }
//            }
//        }
//    }
//}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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

    forAll(surfaceIndex_, facei)
    {
        if (surfaceIndex_[facei] >= 0 && isMasterFace.get(facei) == 1)
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
//    forAll(globalRegion, celli)
//    {
//        Map<label>::const_iterator fndRegion =
//            regionToDist.find(globalRegion[celli]);
//
//        if (fndRegion != regionToDist.end())
//        {
//            distribution[celli] = fndRegion();
//        }
//        else
//        {
//            // region is local to the processor.
//            label localRegionI = globalToLocalRegion[globalRegion[celli]];
//
//            distribution[celli] = regionDistribution[localRegionI];
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
                            label facei = fZone[i];
                            if (blockedFace[facei])
                            {
                                if
                                (
                                    mesh_.isInternalFace(facei)
                                 || pbm[pbm.whichPatch(facei)].coupled()
                                )
                                {
                                    blockedFace[facei] = false;
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
                forAll(separatedCoupledFace, facei)
                {
                    if (separatedCoupledFace[facei])
                    {
                        if (blockedFace[facei])
                        {
                            blockedFace[facei] = false;
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
                forAll(coupledFace, facei)
                {
                    if (coupledFace[facei] != -1 && facei < coupledFace[facei])
                    {
                        couples[nCpl++] = labelPair(facei, coupledFace[facei]);
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

            // if (nUnblocked > 0 || nCouples > 0)
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
            //    forAll(nProcCells, proci)
            //    {
            //        Info<< "    " << proci << '\t' << nProcCells[proci]
            //            << endl;
            //    }
            //    Info<< endl;
            //}
            // else
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
        // else
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
            forAll(nProcCells, proci)
            {
                Pout<< "    " << proci << '\t' << nProcCells[proci] << endl;
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


Foam::labelList Foam::meshRefinement::intersectedFaces() const
{
    label nBoundaryFaces = 0;

    forAll(surfaceIndex_, facei)
    {
        if (surfaceIndex_[facei] != -1)
        {
            nBoundaryFaces++;
        }
    }

    labelList surfaceFaces(nBoundaryFaces);
    nBoundaryFaces = 0;

    forAll(surfaceIndex_, facei)
    {
        if (surfaceIndex_[facei] != -1)
        {
            surfaceFaces[nBoundaryFaces++] = facei;
        }
    }
    return surfaceFaces;
}


Foam::labelList Foam::meshRefinement::intersectedPoints() const
{
    const faceList& faces = mesh_.faces();

    // Mark all points on faces that will become baffles
    PackedBoolList isBoundaryPoint(mesh_.nPoints());
    label nBoundaryPoints = 0;

    forAll(surfaceIndex_, facei)
    {
        if (surfaceIndex_[facei] != -1)
        {
            const face& f = faces[facei];

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
    // labelList adaptPatchIDs(meshedPatches());
    // forAll(adaptPatchIDs, i)
    //{
    //    label patchi = adaptPatchIDs[i];
    //
    //    if (patchi != -1)
    //    {
    //        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
    //
    //        label facei = pp.start();
    //
    //        forAll(pp, i)
    //        {
    //            const face& f = faces[facei];
    //
    //            forAll(f, fp)
    //            {
    //                if (isBoundaryPoint.set(f[fp], 1u))
    //                    nBoundaryPoints++;
    //                }
    //            }
    //            facei++;
    //        }
    //    }
    //}


    // Pack
    labelList boundaryPoints(nBoundaryPoints);
    nBoundaryPoints = 0;
    forAll(isBoundaryPoint, pointi)
    {
        if (isBoundaryPoint.get(pointi) == 1u)
        {
            boundaryPoints[nBoundaryPoints++] = pointi;
        }
    }

    return boundaryPoints;
}


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

        label meshFacei = pp.start();

        forAll(pp, i)
        {
            addressing[nFaces++] = meshFacei++;
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

    forAll(pointPatches, patchi)
    {
        if (isA<processorPointPatch>(pointPatches[patchi]))
        {
            patchFieldTypes[patchi] = calculatedPointPatchVectorField::typeName;
        }
        else if (isA<cyclicPointPatch>(pointPatches[patchi]))
        {
            patchFieldTypes[patchi] = cyclicSlipPointPatchVectorField::typeName;
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
            dimensionedVector("displacement", dimLength, Zero),
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
        forAll(zoneNames, proci)
        {
            if (proci != Pstream::myProcNo())
            {
                if (zoneNames[proci] != zoneNames[Pstream::myProcNo()])
                {
                    FatalErrorInFunction
                        << "faceZones are not synchronised on processors." << nl
                        << "Processor " << proci << " has faceZones "
                        << zoneNames[proci] << nl
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
            label bFacei = fZone[i]-mesh.nInternalFaces();

            if (bFacei >= 0)
            {
                if (faceToZone[bFacei] == -1)
                {
                    faceToZone[bFacei] = zoneI;
                }
                else if (faceToZone[bFacei] == zoneI)
                {
                    FatalErrorInFunction
                        << "Face " << fZone[i] << " in zone "
                        << fZone.name()
                        << " is twice in zone!"
                        << abort(FatalError);
                }
                else
                {
                    FatalErrorInFunction
                        << "Face " << fZone[i] << " in zone "
                        << fZone.name()
                        << " is also in zone "
                        << fZones[faceToZone[bFacei]].name()
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
            FatalErrorInFunction
                << "Face " << mesh.nInternalFaces()+i
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
            small,
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
    forAll(invSumWeight, pointi)
    {
        scalar w = invSumWeight[pointi];

        if (w > 0.0)
        {
            invSumWeight[pointi] = 1.0/w;
        }
    }
}


Foam::label Foam::meshRefinement::appendPatch
(
    fvMesh& mesh,
    const label insertPatchi,
    const word& patchName,
    const dictionary& patchDict
)
{
    // Clear local fields and e.g. polyMesh parallelInfo.
    mesh.clearOut();

    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    label patchi = polyPatches.size();

    // Add polyPatch at the end
    polyPatches.setSize(patchi+1);
    polyPatches.set
    (
        patchi,
        polyPatch::New
        (
            patchName,
            patchDict,
            insertPatchi,
            polyPatches
        )
    );
    fvPatches.setSize(patchi+1);
    fvPatches.set
    (
        patchi,
        fvPatch::New
        (
            polyPatches[patchi],  // point to newly added polyPatch
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
    return patchi;
}


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

    const label patchi = polyPatches.findPatchID(patchName);
    if (patchi != -1)
    {
        // Already there
        return patchi;
    }


    label insertPatchi = polyPatches.size();
    label startFacei = mesh.nFaces();

    forAll(polyPatches, patchi)
    {
        const polyPatch& pp = polyPatches[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            insertPatchi = patchi;
            startFacei = pp.start();
            break;
        }
    }

    dictionary patchDict(patchInfo);
    patchDict.set("nFaces", 0);
    patchDict.set("startFace", startFacei);

    // Below is all quite a hack. Feel free to change once there is a better
    // mechanism to insert and reorder patches.

    label addedPatchi = appendPatch(mesh, insertPatchi, patchName, patchDict);


    // Create reordering list
    // patches before insert position stay as is
    labelList oldToNew(addedPatchi+1);
    for (label i = 0; i < insertPatchi; i++)
    {
        oldToNew[i] = i;
    }
    // patches after insert position move one up
    for (label i = insertPatchi; i < addedPatchi; i++)
    {
        oldToNew[i] = i+1;
    }
    // appended patch gets moved to insert position
    oldToNew[addedPatchi] = insertPatchi;

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

    return insertPatchi;
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
        label patchi = addPatch(mesh_, name, patchInfo);

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
//        label patchi = fvMeshTools::addPatch
//        (
//            mesh_,
//            ppPtr(),
//            dictionary(),       // optional field values
//            calculatedFvPatchField<scalar>::typeName,
//            true
//        );

        // Store
        meshedPatches_.append(name);

        return patchi;
    }
}


Foam::labelList Foam::meshRefinement::meshedPatches() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    DynamicList<label> patchIDs(meshedPatches_.size());
    forAll(meshedPatches_, i)
    {
        label patchi = patches.findPatchID(meshedPatches_[i]);

        if (patchi == -1)
        {
            FatalErrorInFunction
                << "Problem : did not find patch " << meshedPatches_[i]
                << endl << "Valid patches are " << patches.names()
                << abort(FatalError);
        }
        if (!polyPatch::constraintType(patches[patchi].type()))
        {
            patchIDs.append(patchi);
        }
    }

    return patchIDs;
}


void Foam::meshRefinement::selectSeparatedCoupledFaces(boolList& selected) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchi)
    {
        // Check all coupled. Avoid using .coupled() so we also pick up AMI.
        if (isA<coupledPolyPatch>(patches[patchi]))
        {
            const coupledPolyPatch& cpp = refCast<const coupledPolyPatch>
            (
                patches[patchi]
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
    label celli = mesh.findCell(p);
    if (celli != -1)
    {
        regionI = cellToRegion[celli];
    }
    reduce(regionI, maxOp<label>());

    if (regionI == -1)
    {
        // See if we can perturb a bit
        celli = mesh.findCell(p+perturbVec);
        if (celli != -1)
        {
            regionI = cellToRegion[celli];
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
        FatalErrorInFunction
            << "Point " << keepPoint
            << " is not inside the mesh." << nl
            << "Bounding box of the mesh:" << mesh_.bounds()
            << exit(FatalError);
    }

    // Subset
    // ~~~~~~

    // Get cells to remove
    DynamicList<label> cellsToRemove(mesh_.nCells());
    forAll(cellRegion, celli)
    {
        if (cellRegion[celli] != regionI)
        {
            cellsToRemove.append(celli);
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
        // FatalErrorInFunction
        //    << "Removing non-reachable cells should only expose"
        //    << " boundary faces" << nl
        //    << "ExposedFaces:" << exposedFaces << abort(FatalError);

        // Patch for exposed faces for lack of anything sensible.
        label defaultPatch = 0;
        if (globalToMasterPatch.size())
        {
            defaultPatch = globalToMasterPatch[0];
        }

        WarningInFunction
            << "Removing non-reachable cells exposes "
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
        // Get local mesh bounding box. Single box for now.
        List<treeBoundBox> meshBb(1);
        treeBoundBox& bb = meshBb[0];
        bb = treeBoundBox(mesh_.points()).extend(1e-4);

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

            forAll(newFaceData, facei)
            {
                label oldFacei = map.faceMap()[facei];

                if (oldFacei >= 0 && map.reverseFaceMap()[oldFacei] == facei)
                {
                    newFaceData[facei] = data[oldFacei];
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

            forAll(map.faceMap(), facei)
            {
                label oldFacei = map.faceMap()[facei];

                if (oldFacei >= 0)
                {
                    if (reverseFaceMap[oldFacei] != facei)
                    {
                        // facei is slave face. Mark old face.
                        reverseFaceMap[oldFacei] = -1;
                    }
                }
            }

            // 2. Map only faces with intact reverseFaceMap
            labelList newFaceData(map.faceMap().size(), -1);
            forAll(newFaceData, facei)
            {
                label oldFacei = map.faceMap()[facei];

                if (oldFacei >= 0)
                {
                    if (reverseFaceMap[oldFacei] == facei)
                    {
                        newFaceData[facei] = data[oldFacei];
                    }
                }
            }
            data.transfer(newFaceData);
        }
    }
}


bool Foam::meshRefinement::write() const
{
    bool writeOk = mesh_.write();

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
    forAll(meshPoints, pointi)
    {
        myPoints[pointi] = globalPoints.toGlobal(pointi);
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
    forAll(meshPoints, pointi)
    {
        if (myPoints[pointi] == globalPoints.toGlobal(pointi))
        {
            isPatchMasterPoint[pointi] = true;
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


    // if (debug)
    {
        const labelList& cellLevel = meshCutter_.cellLevel();

        labelList nCells(gMax(cellLevel)+1, 0);

        forAll(cellLevel, celli)
        {
            nCells[cellLevel[celli]]++;
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
            dimensionedScalar("zero", dimless, 0)
        );

        const labelList& cellLevel = meshCutter_.cellLevel();

        forAll(volRefLevel, celli)
        {
            volRefLevel[celli] = cellLevel[celli];
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

        forAll(pointRefLevel, pointi)
        {
            pointRefLevel[pointi] = pointLevel[pointi];
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
            label facei = intersectionFaces[i];
            start[i] = cellCentres[mesh_.faceOwner()[facei]];

            if (mesh_.isInternalFace(facei))
            {
                end[i] = cellCentres[mesh_.faceNeighbour()[facei]];
            }
            else
            {
                end[i] = neiCc[facei-mesh_.nInternalFaces()];
            }
        }

        // Extend segments a bit
        {
            const vectorField smallVec(rootSmall*(end-start));
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

    if (writeFlags && !(writeFlags & NOWRITEREFINEMENT))
    {
        meshCutter_.write();
        surfaceIndex_.write();
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
