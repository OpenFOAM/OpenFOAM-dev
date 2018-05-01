/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "fvMeshDistribute.H"
#include "PstreamCombineReduceOps.H"
#include "fvMeshAdder.H"
#include "faceCoupleInfo.H"
#include "processorFvPatchField.H"
#include "processorFvsPatchField.H"
#include "processorCyclicPolyPatch.H"
#include "processorCyclicFvPatchField.H"
#include "polyTopoChange.H"
#include "removeCells.H"
#include "polyModifyFace.H"
#include "polyRemovePoint.H"
#include "mapDistributePolyMesh.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "CompactListList.H"
#include "fvMeshTools.H"
#include "ListOps.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshDistribute, 0);

//- Less function class that can be used for sorting processor patches
class lessProcPatches
{
    const labelList& nbrProc_;
    const labelList& referPatchID_;

public:

    lessProcPatches( const labelList& nbrProc, const labelList& referPatchID)
    :
        nbrProc_(nbrProc),
        referPatchID_(referPatchID)
    {}

    bool operator()(const label a, const label b)
    {
        if (nbrProc_[a] < nbrProc_[b])
        {
            return true;
        }
        else if (nbrProc_[a] > nbrProc_[b])
        {
            return false;
        }
        else
        {
            // Equal neighbour processor
            return referPatchID_[a] < referPatchID_[b];
        }
    }
};

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshDistribute::inplaceRenumberWithFlip
(
    const labelUList& oldToNew,
    const bool oldToNewHasFlip,
    const bool lstHasFlip,
    labelUList& lst
)
{
    if (!lstHasFlip && !oldToNewHasFlip)
    {
        Foam::inplaceRenumber(oldToNew, lst);
    }
    else
    {
        // Either input data or map encodes sign so result encodes sign

        forAll(lst, elemI)
        {
            // Extract old value and sign
            label val = lst[elemI];
            label sign = 1;
            if (lstHasFlip)
            {
                if (val > 0)
                {
                    val = val-1;
                }
                else if (val < 0)
                {
                    val = -val-1;
                    sign = -1;
                }
                else
                {
                    FatalErrorInFunction
                        << "Problem : zero value " << val
                        << " at index " << elemI << " out of " << lst.size()
                        << " list with flip bit" << exit(FatalError);
                }
            }


            // Lookup new value and possibly change sign
            label newVal = oldToNew[val];

            if (oldToNewHasFlip)
            {
                if (newVal > 0)
                {
                    newVal = newVal-1;
                }
                else if (newVal < 0)
                {
                    newVal = -newVal-1;
                    sign = -sign;
                }
                else
                {
                    FatalErrorInFunction
                        << "Problem : zero value " << newVal
                        << " at index " << elemI << " out of "
                        << oldToNew.size()
                        << " list with flip bit" << exit(FatalError);
                }
            }


            // Encode new value and sign
            lst[elemI] = sign*(newVal+1);
        }
    }
}


Foam::labelList Foam::fvMeshDistribute::select
(
    const bool selectEqual,
    const labelList& values,
    const label value
)
{
    label n = 0;

    forAll(values, i)
    {
        if (selectEqual == (values[i] == value))
        {
            n++;
        }
    }

    labelList indices(n);
    n = 0;

    forAll(values, i)
    {
        if (selectEqual == (values[i] == value))
        {
            indices[n++] = i;
        }
    }
    return indices;
}


// Check all procs have same names and in exactly same order.
void Foam::fvMeshDistribute::checkEqualWordList
(
    const string& msg,
    const wordList& lst
)
{
    List<wordList> allNames(Pstream::nProcs());
    allNames[Pstream::myProcNo()] = lst;
    Pstream::gatherList(allNames);
    Pstream::scatterList(allNames);

    for (label proci = 1; proci < Pstream::nProcs(); proci++)
    {
        if (allNames[proci] != allNames[0])
        {
            FatalErrorInFunction
                << "When checking for equal " << msg.c_str() << " :" << endl
                << "processor0 has:" << allNames[0] << endl
                << "processor" << proci << " has:" << allNames[proci] << endl
                << msg.c_str() << " need to be synchronised on all processors."
                << exit(FatalError);
        }
    }
}


Foam::wordList Foam::fvMeshDistribute::mergeWordList(const wordList& procNames)
{
    List<wordList> allNames(Pstream::nProcs());
    allNames[Pstream::myProcNo()] = procNames;
    Pstream::gatherList(allNames);
    Pstream::scatterList(allNames);

    HashSet<word> mergedNames;
    forAll(allNames, proci)
    {
        forAll(allNames[proci], i)
        {
            mergedNames.insert(allNames[proci][i]);
        }
    }
    return mergedNames.toc();
}


// Print some info on mesh.
void Foam::fvMeshDistribute::printMeshInfo(const fvMesh& mesh)
{
    Pout<< "Primitives:" << nl
        << "    points       :" << mesh.nPoints() << nl
        << "    bb           :" << boundBox(mesh.points(), false) << nl
        << "    internalFaces:" << mesh.nInternalFaces() << nl
        << "    faces        :" << mesh.nFaces() << nl
        << "    cells        :" << mesh.nCells() << nl;

    const fvBoundaryMesh& patches = mesh.boundary();

    Pout<< "Patches:" << endl;
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi].patch();

        Pout<< "    " << patchi << " name:" << pp.name()
            << " size:" << pp.size()
            << " start:" << pp.start()
            << " type:" << pp.type()
            << endl;
    }

    if (mesh.pointZones().size())
    {
        Pout<< "PointZones:" << endl;
        forAll(mesh.pointZones(), zoneI)
        {
            const pointZone& pz = mesh.pointZones()[zoneI];
            Pout<< "    " << zoneI << " name:" << pz.name()
                << " size:" << pz.size()
                << endl;
        }
    }
    if (mesh.faceZones().size())
    {
        Pout<< "FaceZones:" << endl;
        forAll(mesh.faceZones(), zoneI)
        {
            const faceZone& fz = mesh.faceZones()[zoneI];
            Pout<< "    " << zoneI << " name:" << fz.name()
                << " size:" << fz.size()
                << endl;
        }
    }
    if (mesh.cellZones().size())
    {
        Pout<< "CellZones:" << endl;
        forAll(mesh.cellZones(), zoneI)
        {
            const cellZone& cz = mesh.cellZones()[zoneI];
            Pout<< "    " << zoneI << " name:" << cz.name()
                << " size:" << cz.size()
                << endl;
        }
    }
}


void Foam::fvMeshDistribute::printCoupleInfo
(
    const primitiveMesh& mesh,
    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourcePatch,
    const labelList& sourceNewNbrProc
)
{
    Pout<< nl
        << "Current coupling info:"
        << endl;

    forAll(sourceFace, bFacei)
    {
        label meshFacei = mesh.nInternalFaces() + bFacei;

        Pout<< "    meshFace:" << meshFacei
            << " fc:" << mesh.faceCentres()[meshFacei]
            << " connects to proc:" << sourceProc[bFacei]
            << "/face:" << sourceFace[bFacei]
            << " which will move to proc:" << sourceNewNbrProc[bFacei]
            << endl;
    }
}


// Finds (non-empty) patch that exposed internal and proc faces can be put into.
Foam::label Foam::fvMeshDistribute::findNonEmptyPatch() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label nonEmptyPatchi = -1;

    forAllReverse(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (!isA<emptyPolyPatch>(pp) && !pp.coupled())
        {
            nonEmptyPatchi = patchi;
            break;
        }
    }

    if (nonEmptyPatchi == -1)
    {
        FatalErrorInFunction
            << "Cannot find a patch which is neither of type empty nor"
            << " coupled in patches " << patches.names() << endl
            << "There has to be at least one such patch for"
            << " distribution to work" << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "findNonEmptyPatch : using patch " << nonEmptyPatchi
            << " name:" << patches[nonEmptyPatchi].name()
            << " type:" << patches[nonEmptyPatchi].type()
            << " to put exposed faces into." << endl;
    }


    // Do additional test for processor patches intermingled with non-proc
    // patches.
    label procPatchi = -1;

    forAll(patches, patchi)
    {
        if (isA<processorPolyPatch>(patches[patchi]))
        {
            procPatchi = patchi;
        }
        else if (procPatchi != -1)
        {
            FatalErrorInFunction
                << "Processor patches should be at end of patch list."
                << endl
                << "Have processor patch " << procPatchi
                << " followed by non-processor patch " << patchi
                << " in patches " << patches.names()
                << abort(FatalError);
        }
    }

    return nonEmptyPatchi;
}


// Delete all processor patches. Move any processor faces into the last
// non-processor patch.
Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::deleteProcPatches
(
    const label destinationPatch
)
{
    // New patchID per boundary faces to be repatched. Is -1 (no change)
    // or new patchID
    labelList newPatchID(mesh_.nFaces() - mesh_.nInternalFaces(), -1);

    label nProcPatches = 0;

    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            if (debug)
            {
                Pout<< "Moving all faces of patch " << pp.name()
                    << " into patch " << destinationPatch
                    << endl;
            }

            label offset = pp.start() - mesh_.nInternalFaces();

            forAll(pp, i)
            {
                newPatchID[offset+i] = destinationPatch;
            }

            nProcPatches++;
        }
    }

    // Note: order of boundary faces been kept the same since the
    // destinationPatch is at the end and we have visited the patches in
    // incremental order.
    labelListList dummyFaceMaps;
    autoPtr<mapPolyMesh> map = repatch(newPatchID, dummyFaceMaps);


    // Delete (now empty) processor patches.
    {
        labelList oldToNew(identity(mesh_.boundaryMesh().size()));
        label newI = 0;
        // Non processor patches first
        forAll(mesh_.boundaryMesh(), patchi)
        {
            if (!isA<processorPolyPatch>(mesh_.boundaryMesh()[patchi]))
            {
                oldToNew[patchi] = newI++;
            }
        }
        label nNonProcPatches = newI;

        // Processor patches as last
        forAll(mesh_.boundaryMesh(), patchi)
        {
            if (isA<processorPolyPatch>(mesh_.boundaryMesh()[patchi]))
            {
                oldToNew[patchi] = newI++;
            }
        }
        fvMeshTools::reorderPatches(mesh_, oldToNew, nNonProcPatches, false);
    }

    return map;
}


// Repatch the mesh.
Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::repatch
(
    const labelList& newPatchID,         // per boundary face -1 or new patchID
    labelListList& constructFaceMap
)
{
    polyTopoChange meshMod(mesh_);

    forAll(newPatchID, bFacei)
    {
        if (newPatchID[bFacei] != -1)
        {
            label facei = mesh_.nInternalFaces() + bFacei;

            label zoneID = mesh_.faceZones().whichZone(facei);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = mesh_.faceZones()[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh_.faces()[facei],       // modified face
                    facei,                      // label of face
                    mesh_.faceOwner()[facei],   // owner
                    -1,                         // neighbour
                    false,                      // face flip
                    newPatchID[bFacei],         // patch for face
                    false,                      // remove from zone
                    zoneID,                     // zone for face
                    zoneFlip                    // face flip in zone
                )
            );
        }
    }


    // Do mapping of fields from one patchField to the other ourselves since
    // is currently not supported by updateMesh.

    // Store boundary fields (we only do this for surfaceFields)
    PtrList<FieldField<fvsPatchField, scalar>> sFlds;
    saveBoundaryFields<scalar, surfaceMesh>(sFlds);
    PtrList<FieldField<fvsPatchField, vector>> vFlds;
    saveBoundaryFields<vector, surfaceMesh>(vFlds);
    PtrList<FieldField<fvsPatchField, sphericalTensor>> sptFlds;
    saveBoundaryFields<sphericalTensor, surfaceMesh>(sptFlds);
    PtrList<FieldField<fvsPatchField, symmTensor>> sytFlds;
    saveBoundaryFields<symmTensor, surfaceMesh>(sytFlds);
    PtrList<FieldField<fvsPatchField, tensor>> tFlds;
    saveBoundaryFields<tensor, surfaceMesh>(tFlds);

    // Change the mesh (no inflation). Note: parallel comms allowed.
    //
    // NOTE: there is one very particular problem with this ordering.
    // We first create the processor patches and use these to merge out
    // shared points (see mergeSharedPoints below). So temporarily points
    // and edges do not match!

    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields. No inflation, parallel sync.
    mesh_.updateMesh(map);

    // Map patch fields using stored boundary fields. Note: assumes order
    // of fields has not changed in object registry!
    mapBoundaryFields<scalar, surfaceMesh>(map, sFlds);
    mapBoundaryFields<vector, surfaceMesh>(map, vFlds);
    mapBoundaryFields<sphericalTensor, surfaceMesh>(map, sptFlds);
    mapBoundaryFields<symmTensor, surfaceMesh>(map, sytFlds);
    mapBoundaryFields<tensor, surfaceMesh>(map, tFlds);


    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }

    // Adapt constructMaps.

    if (debug)
    {
        label index = findIndex(map().reverseFaceMap(), -1);

        if (index != -1)
        {
            FatalErrorInFunction
                << "reverseFaceMap contains -1 at index:"
                << index << endl
                << "This means that the repatch operation was not just"
                << " a shuffle?" << abort(FatalError);
        }
    }

    forAll(constructFaceMap, proci)
    {
        inplaceRenumberWithFlip
        (
            map().reverseFaceMap(),
            false,
            true,
            constructFaceMap[proci]
        );
    }


    return map;
}


// Detect shared points. Need processor patches to be present.
// Background: when adding bits of mesh one can get points which
// share the same position but are only detectable to be topologically
// the same point when doing parallel analysis. This routine will
// merge those points.
Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::mergeSharedPoints
(
    const labelList& pointToGlobalMaster,
    labelListList& constructPointMap
)
{
    // Find out which sets of points get merged and create a map from
    // mesh point to unique point.

    label nShared = 0;
    forAll(pointToGlobalMaster, pointi)
    {
        if (pointToGlobalMaster[pointi] != -1)
        {
            nShared++;
        }
    }

    Map<label> globalMasterToLocalMaster(2*nShared);
    Map<label> pointToMaster(2*nShared);

    forAll(pointToGlobalMaster, pointi)
    {
        label globali = pointToGlobalMaster[pointi];
        if (globali != -1)
        {
            Map<label>::const_iterator iter = globalMasterToLocalMaster.find
            (
                globali
            );

            if (iter == globalMasterToLocalMaster.end())
            {
                // Found first point. Designate as master
                globalMasterToLocalMaster.insert(globali, pointi);
                pointToMaster.insert(pointi, pointi);
            }
            else
            {
                pointToMaster.insert(pointi, iter());
            }
        }
    }

    if (returnReduce(pointToMaster.size(), sumOp<label>()) == 0)
    {
        return autoPtr<mapPolyMesh>(nullptr);
    }


    polyTopoChange meshMod(mesh_);

    fvMeshAdder::mergePoints(mesh_, pointToMaster, meshMod);

    // Change the mesh (no inflation). Note: parallel comms allowed.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields. No inflation, parallel sync.
    mesh_.updateMesh(map);

    // Adapt constructMaps for merged points.
    forAll(constructPointMap, proci)
    {
        labelList& constructMap = constructPointMap[proci];

        forAll(constructMap, i)
        {
            label oldPointi = constructMap[i];

            label newPointi = map().reversePointMap()[oldPointi];

            if (newPointi < -1)
            {
                constructMap[i] = -newPointi-2;
            }
            else if (newPointi >= 0)
            {
                constructMap[i] = newPointi;
            }
            else
            {
                FatalErrorInFunction
                    << "Problem. oldPointi:" << oldPointi
                    << " newPointi:" << newPointi << abort(FatalError);
            }
        }
    }
    return map;
}


void Foam::fvMeshDistribute::getCouplingData
(
    const labelList& distribution,
    labelList& sourceFace,
    labelList& sourceProc,
    labelList& sourcePatch,
    labelList& sourceNewNbrProc,
    labelList& sourcePointMaster
) const
{
    // Construct the coupling information for all (boundary) faces and
    // points

    label nBnd = mesh_.nFaces() - mesh_.nInternalFaces();
    sourceFace.setSize(nBnd);
    sourceProc.setSize(nBnd);
    sourcePatch.setSize(nBnd);
    sourceNewNbrProc.setSize(nBnd);

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Get neighbouring meshFace labels and new processor of coupled boundaries.
    labelList nbrFaces(nBnd, -1);
    labelList nbrNewNbrProc(nBnd, -1);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            label offset = pp.start() - mesh_.nInternalFaces();

            // Mesh labels of faces on this side
            forAll(pp, i)
            {
                label bndI = offset + i;
                nbrFaces[bndI] = pp.start()+i;
            }

            // Which processor they will end up on
            SubList<label>(nbrNewNbrProc, pp.size(), offset) =
                UIndirectList<label>(distribution, pp.faceCells())();
        }
    }


    // Exchange the boundary data
    syncTools::swapBoundaryFaceList(mesh_, nbrFaces);
    syncTools::swapBoundaryFaceList(mesh_, nbrNewNbrProc);


    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        label offset = pp.start() - mesh_.nInternalFaces();

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            // Check which of the two faces we store.

            if (procPatch.owner())
            {
                // Use my local face labels
                forAll(pp, i)
                {
                    label bndI = offset + i;
                    sourceFace[bndI] = pp.start()+i;
                    sourceProc[bndI] = Pstream::myProcNo();
                    sourceNewNbrProc[bndI] = nbrNewNbrProc[bndI];
                }
            }
            else
            {
                // Use my neighbours face labels
                forAll(pp, i)
                {
                    label bndI = offset + i;
                    sourceFace[bndI] = nbrFaces[bndI];
                    sourceProc[bndI] = procPatch.neighbProcNo();
                    sourceNewNbrProc[bndI] = nbrNewNbrProc[bndI];
                }
            }


            label patchi = -1;
            if (isA<processorCyclicPolyPatch>(pp))
            {
                patchi = refCast<const processorCyclicPolyPatch>
                (
                    pp
                ).referPatchID();
            }

            forAll(pp, i)
            {
                label bndI = offset + i;
                sourcePatch[bndI] = patchi;
            }
        }
        else if (isA<cyclicPolyPatch>(pp))
        {
            const cyclicPolyPatch& cpp = refCast<const cyclicPolyPatch>(pp);

            if (cpp.owner())
            {
                forAll(pp, i)
                {
                    label bndI = offset + i;
                    sourceFace[bndI] = pp.start()+i;
                    sourceProc[bndI] = Pstream::myProcNo();
                    sourcePatch[bndI] = patchi;
                    sourceNewNbrProc[bndI] = nbrNewNbrProc[bndI];
                }
            }
            else
            {
                forAll(pp, i)
                {
                    label bndI = offset + i;
                    sourceFace[bndI] = nbrFaces[bndI];
                    sourceProc[bndI] = Pstream::myProcNo();
                    sourcePatch[bndI] = patchi;
                    sourceNewNbrProc[bndI] = nbrNewNbrProc[bndI];
                }
            }
        }
        else
        {
            // Normal (physical) boundary
            forAll(pp, i)
            {
                label bndI = offset + i;
                sourceFace[bndI] = -1;
                sourceProc[bndI] = -1;
                sourcePatch[bndI] = patchi;
                sourceNewNbrProc[bndI] = -1;
            }
        }
    }


    // Collect coupled (collocated) points
    sourcePointMaster.setSize(mesh_.nPoints());
    sourcePointMaster = -1;
    {
        // Assign global master point
        const globalIndex globalPoints(mesh_.nPoints());

        const globalMeshData& gmd = mesh_.globalData();
        const indirectPrimitivePatch& cpp = gmd.coupledPatch();
        const labelList& meshPoints = cpp.meshPoints();
        const mapDistribute& slavesMap = gmd.globalCoPointSlavesMap();
        const labelListList& slaves = gmd.globalCoPointSlaves();

        labelList elems(slavesMap.constructSize(), -1);
        forAll(meshPoints, pointi)
        {
            const labelList& slots = slaves[pointi];

            if (slots.size())
            {
                // pointi is a master. Assign a unique label.

                label globalPointi = globalPoints.toGlobal(meshPoints[pointi]);
                elems[pointi] = globalPointi;
                forAll(slots, i)
                {
                    label sloti = slots[i];
                    if (sloti >= meshPoints.size())
                    {
                        // Filter out local collocated points. We don't want
                        // to merge these
                        elems[slots[i]] = globalPointi;
                    }
                }
            }
        }

        // Push slave-slot data back to slaves
        slavesMap.reverseDistribute(elems.size(), elems, false);

        // Extract back onto mesh
        forAll(meshPoints, pointi)
        {
            sourcePointMaster[meshPoints[pointi]] = elems[pointi];
        }
    }
}


// Subset the neighbourCell/neighbourProc fields
void Foam::fvMeshDistribute::subsetCouplingData
(
    const fvMesh& mesh,
    const labelList& pointMap,
    const labelList& faceMap,
    const labelList& cellMap,

    const labelList& oldDistribution,
    const labelList& oldFaceOwner,
    const labelList& oldFaceNeighbour,
    const label oldInternalFaces,

    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourcePatch,
    const labelList& sourceNewNbrProc,
    const labelList& sourcePointMaster,

    labelList& subFace,
    labelList& subProc,
    labelList& subPatch,
    labelList& subNewNbrProc,
    labelList& subPointMaster
)
{
    subFace.setSize(mesh.nFaces() - mesh.nInternalFaces());
    subProc.setSize(mesh.nFaces() - mesh.nInternalFaces());
    subPatch.setSize(mesh.nFaces() - mesh.nInternalFaces());
    subNewNbrProc.setSize(mesh.nFaces() - mesh.nInternalFaces());

    forAll(subFace, newBFacei)
    {
        label newFacei = newBFacei + mesh.nInternalFaces();

        label oldFacei = faceMap[newFacei];

        // Was oldFacei internal face? If so which side did we get.
        if (oldFacei < oldInternalFaces)
        {
            subFace[newBFacei] = oldFacei;
            subProc[newBFacei] = Pstream::myProcNo();
            subPatch[newBFacei] = -1;

            label oldOwn = oldFaceOwner[oldFacei];
            label oldNei = oldFaceNeighbour[oldFacei];

            if (oldOwn == cellMap[mesh.faceOwner()[newFacei]])
            {
                // We kept the owner side. Where does the neighbour move to?
                subNewNbrProc[newBFacei] = oldDistribution[oldNei];
            }
            else
            {
                // We kept the neighbour side.
                subNewNbrProc[newBFacei] = oldDistribution[oldOwn];
            }
        }
        else
        {
            // Was boundary face. Take over boundary information
            label oldBFacei = oldFacei - oldInternalFaces;

            subFace[newBFacei] = sourceFace[oldBFacei];
            subProc[newBFacei] = sourceProc[oldBFacei];
            subPatch[newBFacei] = sourcePatch[oldBFacei];
            subNewNbrProc[newBFacei] = sourceNewNbrProc[oldBFacei];
        }
    }


    subPointMaster = UIndirectList<label>(sourcePointMaster, pointMap);
}


// Find cells on mesh whose faceID/procID match the neighbour cell/proc of
// domainMesh. Store the matching face.
void Foam::fvMeshDistribute::findCouples
(
    const primitiveMesh& mesh,
    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourcePatch,

    const label domain,
    const primitiveMesh& domainMesh,
    const labelList& domainFace,
    const labelList& domainProc,
    const labelList& domainPatch,

    labelList& masterCoupledFaces,
    labelList& slaveCoupledFaces
)
{
    // Store domain neighbour as map so we can easily look for pair
    // with same face+proc.
    HashTable<label, labelPair, labelPair::Hash<>> map(domainFace.size());

    forAll(domainProc, bFacei)
    {
        if (domainProc[bFacei] != -1 && domainPatch[bFacei] == -1)
        {
            map.insert
            (
                labelPair(domainFace[bFacei], domainProc[bFacei]),
                bFacei
            );
        }
    }


    // Try to match mesh data.

    masterCoupledFaces.setSize(domainFace.size());
    slaveCoupledFaces.setSize(domainFace.size());
    label coupledI = 0;

    forAll(sourceFace, bFacei)
    {
        if (sourceProc[bFacei] != -1 && sourcePatch[bFacei] == -1)
        {
            labelPair myData(sourceFace[bFacei], sourceProc[bFacei]);

            HashTable<label, labelPair, labelPair::Hash<>>::const_iterator
                iter = map.find(myData);

            if (iter != map.end())
            {
                label nbrBFacei = iter();

                masterCoupledFaces[coupledI] = mesh.nInternalFaces() + bFacei;
                slaveCoupledFaces[coupledI] =
                    domainMesh.nInternalFaces()
                  + nbrBFacei;

                coupledI++;
            }
        }
    }

    masterCoupledFaces.setSize(coupledI);
    slaveCoupledFaces.setSize(coupledI);

    if (debug)
    {
        Pout<< "findCouples : found " << coupledI
            << " faces that will be stitched" << nl << endl;
    }
}


// Map data on boundary faces to new mesh (resulting from adding two meshes)
Foam::labelList Foam::fvMeshDistribute::mapBoundaryData
(
    const primitiveMesh& mesh,      // mesh after adding
    const mapAddedPolyMesh& map,
    const labelList& boundaryData0, // on mesh before adding
    const label nInternalFaces1,
    const labelList& boundaryData1  // on added mesh
)
{
    labelList newBoundaryData(mesh.nFaces() - mesh.nInternalFaces());

    forAll(boundaryData0, oldBFacei)
    {
        label newFacei = map.oldFaceMap()[oldBFacei + map.nOldInternalFaces()];

        // Face still exists (is necessary?) and still boundary face
        if (newFacei >= 0 && newFacei >= mesh.nInternalFaces())
        {
            newBoundaryData[newFacei - mesh.nInternalFaces()] =
                boundaryData0[oldBFacei];
        }
    }

    forAll(boundaryData1, addedBFacei)
    {
        label newFacei = map.addedFaceMap()[addedBFacei + nInternalFaces1];

        if (newFacei >= 0 && newFacei >= mesh.nInternalFaces())
        {
            newBoundaryData[newFacei - mesh.nInternalFaces()] =
                boundaryData1[addedBFacei];
        }
    }

    return newBoundaryData;
}


Foam::labelList Foam::fvMeshDistribute::mapPointData
(
    const primitiveMesh& mesh,      // mesh after adding
    const mapAddedPolyMesh& map,
    const labelList& boundaryData0, // on mesh before adding
    const labelList& boundaryData1  // on added mesh
)
{
    labelList newBoundaryData(mesh.nPoints());

    forAll(boundaryData0, oldPointi)
    {
        label newPointi = map.oldPointMap()[oldPointi];

        // Point still exists (is necessary?)
        if (newPointi >= 0)
        {
            newBoundaryData[newPointi] = boundaryData0[oldPointi];
        }
    }

    forAll(boundaryData1, addedPointi)
    {
        label newPointi = map.addedPointMap()[addedPointi];

        if (newPointi >= 0)
        {
            newBoundaryData[newPointi] = boundaryData1[addedPointi];
        }
    }

    return newBoundaryData;
}


// Remove cells. Add all exposed faces to patch oldInternalPatchi
Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::doRemoveCells
(
    const labelList& cellsToRemove,
    const label oldInternalPatchi
)
{
    // Mesh change engine
    polyTopoChange meshMod(mesh_);

    // Cell removal topo engine. Do NOT synchronize parallel since
    // we are doing a local cell removal.
    removeCells cellRemover(mesh_, false);

    // Get all exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));

    // Insert the topo changes
    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        labelList(exposedFaces.size(), oldInternalPatchi),  // patch for exposed
                                                            // faces.
        meshMod
    );


    //// Generate test field
    // tmp<surfaceScalarField> sfld(generateTestField(mesh_));

    // Save internal fields (note: not as DimensionedFields since would
    // get mapped)
    PtrList<Field<scalar>> sFlds;
    saveInternalFields(sFlds);
    PtrList<Field<vector>> vFlds;
    saveInternalFields(vFlds);
    PtrList<Field<sphericalTensor>> sptFlds;
    saveInternalFields(sptFlds);
    PtrList<Field<symmTensor>> sytFlds;
    saveInternalFields(sytFlds);
    PtrList<Field<tensor>> tFlds;
    saveInternalFields(tFlds);

    // Change the mesh. No inflation. Note: no parallel comms allowed.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, false);

    // Update fields
    mesh_.updateMesh(map);


    // Any exposed faces in a surfaceField will not be mapped. Map the value
    // of these separately (until there is support in all PatchFields for
    // mapping from internal faces ...)

    mapExposedFaces(map(), sFlds);
    mapExposedFaces(map(), vFlds);
    mapExposedFaces(map(), sptFlds);
    mapExposedFaces(map(), sytFlds);
    mapExposedFaces(map(), tFlds);


    //// Test test field
    // testField(sfld);


    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }

    return map;
}


// Delete and add processor patches. Changes mesh and returns per neighbour proc
// the processor patchID.
void Foam::fvMeshDistribute::addProcPatches
(
    const labelList& nbrProc,       // processor that neighbour is now on
    const labelList& referPatchID,  // patchID (or -1) I originated from
    List<Map<label>>& procPatchID
)
{
    // Now use the neighbourFace/Proc to repatch the mesh. These lists
    // contain for all current boundary faces the global patchID (for non-proc
    // patch) or the processor.

    // Determine a visit order such that the processor patches get added
    // in order of increasing neighbour processor (and for same neighbour
    // processor (in case of processor cyclics) in order of increasing
    // 'refer' patch)
    labelList indices;
    sortedOrder(nbrProc, indices, lessProcPatches(nbrProc, referPatchID));

    procPatchID.setSize(Pstream::nProcs());

    forAll(indices, i)
    {
        label bFacei = indices[i];
        label proci = nbrProc[bFacei];

        if (proci != -1 && proci != Pstream::myProcNo())
        {
            if (!procPatchID[proci].found(referPatchID[bFacei]))
            {
                // No patch for neighbour yet. Is either a normal processor
                // patch or a processorCyclic patch.

                if (referPatchID[bFacei] == -1)
                {
                    // Ordinary processor boundary

                    processorPolyPatch pp
                    (
                        0,              // size
                        mesh_.nFaces(),
                        mesh_.boundaryMesh().size(),
                        mesh_.boundaryMesh(),
                        Pstream::myProcNo(),
                        proci
                    );

                    procPatchID[proci].insert
                    (
                        referPatchID[bFacei],
                        fvMeshTools::addPatch
                        (
                            mesh_,
                            pp,
                            dictionary(),   // optional per field patchField
                            processorFvPatchField<scalar>::typeName,
                            false           // not parallel sync
                        )
                    );
                }
                else
                {
                    const coupledPolyPatch& pcPatch
                        = refCast<const coupledPolyPatch>
                          (
                              mesh_.boundaryMesh()[referPatchID[bFacei]]
                          );
                    processorCyclicPolyPatch pp
                    (
                        0,              // size
                        mesh_.nFaces(),
                        mesh_.boundaryMesh().size(),
                        mesh_.boundaryMesh(),
                        Pstream::myProcNo(),
                        proci,
                        pcPatch.name(),
                        pcPatch.transform()
                    );

                    procPatchID[proci].insert
                    (
                        referPatchID[bFacei],
                        fvMeshTools::addPatch
                        (
                            mesh_,
                            pp,
                            dictionary(),   // optional per field patchField
                            processorCyclicFvPatchField<scalar>::typeName,
                            false           // not parallel sync
                        )
                    );
                }
            }
        }
    }
}


// Get boundary faces to be repatched. Is -1 or new patchID
Foam::labelList Foam::fvMeshDistribute::getBoundaryPatch
(
    const labelList& nbrProc,               // new processor per boundary face
    const labelList& referPatchID,          // patchID (or -1) I originated from
    const List<Map<label>>& procPatchID    // per proc the new procPatches
)
{
    labelList patchIDs(nbrProc);

    forAll(nbrProc, bFacei)
    {
        if (nbrProc[bFacei] == Pstream::myProcNo())
        {
            label origPatchi = referPatchID[bFacei];
            patchIDs[bFacei] = origPatchi;
        }
        else if (nbrProc[bFacei] != -1)
        {
            label origPatchi = referPatchID[bFacei];
            patchIDs[bFacei] = procPatchID[nbrProc[bFacei]][origPatchi];
        }
        else
        {
            patchIDs[bFacei] = -1;
        }
    }
    return patchIDs;
}


// Send mesh and coupling data.
void Foam::fvMeshDistribute::sendMesh
(
    const label domain,
    const fvMesh& mesh,

    const wordList& pointZoneNames,
    const wordList& faceZoneNames,
    const wordList& cellZoneNames,

    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourcePatch,
    const labelList& sourceNewNbrProc,
    const labelList& sourcePointMaster,
    Ostream& toDomain
)
{
    if (debug)
    {
        Pout<< "Sending to domain " << domain << nl
            << "    nPoints:" << mesh.nPoints() << nl
            << "    nFaces:" << mesh.nFaces() << nl
            << "    nCells:" << mesh.nCells() << nl
            << "    nPatches:" << mesh.boundaryMesh().size() << nl
            << endl;
    }

    // Assume sparse, possibly overlapping point zones. Get contents
    // in merged-zone indices.
    CompactListList<label> zonePoints;
    {
        const pointZoneMesh& pointZones = mesh.pointZones();

        labelList rowSizes(pointZoneNames.size(), 0);

        forAll(pointZoneNames, nameI)
        {
            label myZoneID = pointZones.findZoneID(pointZoneNames[nameI]);

            if (myZoneID != -1)
            {
                rowSizes[nameI] = pointZones[myZoneID].size();
            }
        }
        zonePoints.setSize(rowSizes);

        forAll(pointZoneNames, nameI)
        {
            label myZoneID = pointZones.findZoneID(pointZoneNames[nameI]);

            if (myZoneID != -1)
            {
                zonePoints[nameI].deepCopy(pointZones[myZoneID]);
            }
        }
    }

    // Assume sparse, possibly overlapping face zones
    CompactListList<label> zoneFaces;
    CompactListList<bool> zoneFaceFlip;
    {
        const faceZoneMesh& faceZones = mesh.faceZones();

        labelList rowSizes(faceZoneNames.size(), 0);

        forAll(faceZoneNames, nameI)
        {
            label myZoneID = faceZones.findZoneID(faceZoneNames[nameI]);

            if (myZoneID != -1)
            {
                rowSizes[nameI] = faceZones[myZoneID].size();
            }
        }

        zoneFaces.setSize(rowSizes);
        zoneFaceFlip.setSize(rowSizes);

        forAll(faceZoneNames, nameI)
        {
            label myZoneID = faceZones.findZoneID(faceZoneNames[nameI]);

            if (myZoneID != -1)
            {
                zoneFaces[nameI].deepCopy(faceZones[myZoneID]);
                zoneFaceFlip[nameI].deepCopy(faceZones[myZoneID].flipMap());
            }
        }
    }

    // Assume sparse, possibly overlapping cell zones
    CompactListList<label> zoneCells;
    {
        const cellZoneMesh& cellZones = mesh.cellZones();

        labelList rowSizes(cellZoneNames.size(), 0);

        forAll(cellZoneNames, nameI)
        {
            label myZoneID = cellZones.findZoneID(cellZoneNames[nameI]);

            if (myZoneID != -1)
            {
                rowSizes[nameI] = cellZones[myZoneID].size();
            }
        }

        zoneCells.setSize(rowSizes);

        forAll(cellZoneNames, nameI)
        {
            label myZoneID = cellZones.findZoneID(cellZoneNames[nameI]);

            if (myZoneID != -1)
            {
                zoneCells[nameI].deepCopy(cellZones[myZoneID]);
            }
        }
    }
    ////- Assume full cell zones
    // labelList cellZoneID;
    // if (hasCellZones)
    //{
    //    cellZoneID.setSize(mesh.nCells());
    //    cellZoneID = -1;
    //
    //    const cellZoneMesh& cellZones = mesh.cellZones();
    //
    //    forAll(cellZones, zoneI)
    //    {
    //        UIndirectList<label>(cellZoneID, cellZones[zoneI]) = zoneI;
    //    }
    //}

    // Send
    toDomain
        << mesh.points()
        << CompactListList<label, face>(mesh.faces())
        << mesh.faceOwner()
        << mesh.faceNeighbour()
        << mesh.boundaryMesh()

        << zonePoints
        << zoneFaces
        << zoneFaceFlip
        << zoneCells

        << sourceFace
        << sourceProc
        << sourcePatch
        << sourceNewNbrProc
        << sourcePointMaster;


    if (debug)
    {
        Pout<< "Started sending mesh to domain " << domain
            << endl;
    }
}


// Receive mesh. Opposite of sendMesh
Foam::autoPtr<Foam::fvMesh> Foam::fvMeshDistribute::receiveMesh
(
    const label domain,
    const wordList& pointZoneNames,
    const wordList& faceZoneNames,
    const wordList& cellZoneNames,
    const Time& runTime,
    labelList& domainSourceFace,
    labelList& domainSourceProc,
    labelList& domainSourcePatch,
    labelList& domainSourceNewNbrProc,
    labelList& domainSourcePointMaster,
    Istream& fromNbr
)
{
    pointField domainPoints(fromNbr);
    faceList domainFaces = CompactListList<label, face>(fromNbr)();
    labelList domainAllOwner(fromNbr);
    labelList domainAllNeighbour(fromNbr);
    PtrList<entry> patchEntries(fromNbr);

    CompactListList<label> zonePoints(fromNbr);
    CompactListList<label> zoneFaces(fromNbr);
    CompactListList<bool> zoneFaceFlip(fromNbr);
    CompactListList<label> zoneCells(fromNbr);

    fromNbr
        >> domainSourceFace
        >> domainSourceProc
        >> domainSourcePatch
        >> domainSourceNewNbrProc
        >> domainSourcePointMaster;

    // Construct fvMesh
    autoPtr<fvMesh> domainMeshPtr
    (
        new fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::NO_READ
            ),
            xferMove(domainPoints),
            xferMove(domainFaces),
            xferMove(domainAllOwner),
            xferMove(domainAllNeighbour),
            false                   // no parallel comms
        )
    );
    fvMesh& domainMesh = domainMeshPtr();

    List<polyPatch*> patches(patchEntries.size());

    forAll(patchEntries, patchi)
    {
        patches[patchi] = polyPatch::New
        (
            patchEntries[patchi].keyword(),
            patchEntries[patchi].dict(),
            patchi,
            domainMesh.boundaryMesh()
        ).ptr();
    }
    // Add patches; no parallel comms
    domainMesh.addFvPatches(patches, false);

    // Construct zones
    List<pointZone*> pZonePtrs(pointZoneNames.size());
    forAll(pZonePtrs, i)
    {
        pZonePtrs[i] = new pointZone
        (
            pointZoneNames[i],
            zonePoints[i],
            i,
            domainMesh.pointZones()
        );
    }

    List<faceZone*> fZonePtrs(faceZoneNames.size());
    forAll(fZonePtrs, i)
    {
        fZonePtrs[i] = new faceZone
        (
            faceZoneNames[i],
            zoneFaces[i],
            zoneFaceFlip[i],
            i,
            domainMesh.faceZones()
        );
    }

    List<cellZone*> cZonePtrs(cellZoneNames.size());
    forAll(cZonePtrs, i)
    {
        cZonePtrs[i] = new cellZone
        (
            cellZoneNames[i],
            zoneCells[i],
            i,
            domainMesh.cellZones()
        );
    }
    domainMesh.addZones(pZonePtrs, fZonePtrs, cZonePtrs);

    return domainMeshPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::fvMeshDistribute::fvMeshDistribute(fvMesh& mesh, const scalar mergeTol)
:
    mesh_(mesh),
    mergeTol_(mergeTol)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::fvMeshDistribute::countCells
(
    const labelList& distribution
)
{
    labelList nCells(Pstream::nProcs(), 0);
    forAll(distribution, celli)
    {
        label newProc = distribution[celli];

        if (newProc < 0 || newProc >= Pstream::nProcs())
        {
            FatalErrorInFunction
                << "Distribution should be in range 0.." << Pstream::nProcs()-1
                << endl
                << "At index " << celli << " distribution:" << newProc
                << abort(FatalError);
        }
        nCells[newProc]++;
    }
    return nCells;
}


Foam::autoPtr<Foam::mapDistributePolyMesh> Foam::fvMeshDistribute::distribute
(
    const labelList& distribution
)
{
    // Some checks on distribution
    if (distribution.size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "Size of distribution:"
            << distribution.size() << " mesh nCells:" << mesh_.nCells()
            << abort(FatalError);
    }


    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Check all processors have same non-proc patches in same order.
    if (patches.checkParallelSync(true))
    {
        FatalErrorInFunction
            << "This application requires all non-processor patches"
            << " to be present in the same order on all patches" << nl
            << "followed by the processor patches (which of course are unique)."
            << nl
            << "Local patches:" << mesh_.boundaryMesh().names()
            << abort(FatalError);
    }

    // Save some data for mapping later on
    const label nOldPoints(mesh_.nPoints());
    const label nOldFaces(mesh_.nFaces());
    const label nOldCells(mesh_.nCells());
    labelList oldPatchStarts(patches.size());
    labelList oldPatchNMeshPoints(patches.size());
    forAll(patches, patchi)
    {
        oldPatchStarts[patchi] = patches[patchi].start();
        oldPatchNMeshPoints[patchi] = patches[patchi].nPoints();
    }


    // Short circuit trivial case.
    if (!Pstream::parRun())
    {
        // Collect all maps and return
        return autoPtr<mapDistributePolyMesh>
        (
            new mapDistributePolyMesh
            (
                mesh_,

                nOldPoints,
                nOldFaces,
                nOldCells,
                oldPatchStarts.xfer(),
                oldPatchNMeshPoints.xfer(),

                labelListList(1, identity(mesh_.nPoints())).xfer(),
                labelListList(1, identity(mesh_.nFaces())).xfer(),
                labelListList(1, identity(mesh_.nCells())).xfer(),
                labelListList(1, identity(patches.size())).xfer(),

                labelListList(1, identity(mesh_.nPoints())).xfer(),
                labelListList(1, identity(mesh_.nFaces())).xfer(),
                labelListList(1, identity(mesh_.nCells())).xfer(),
                labelListList(1, identity(patches.size())).xfer()
            )
        );
    }


    // Collect any zone names
    const wordList pointZoneNames(mergeWordList(mesh_.pointZones().names()));
    const wordList faceZoneNames(mergeWordList(mesh_.faceZones().names()));
    const wordList cellZoneNames(mergeWordList(mesh_.cellZones().names()));


    // Local environment of all boundary faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // A face is uniquely defined by
    //  - proc
    //  - local face no
    //
    // To glue the parts of meshes which can get sent from anywhere we
    // need to know on boundary faces what the above tuple on both sides is.
    // So we need to maintain:
    //  - original face
    //  - original processor id (= trivial)
    // For coupled boundaries (where the faces are 'duplicate') we take the
    // lowest numbered processor as the data to store.
    //
    // Additionally to create the procboundaries we need to know where the owner
    // cell on the other side now is: newNeighbourProc.
    //

    // physical boundary:
    //     sourceProc = -1
    //     sourceNewNbrProc = -1
    //     sourceFace = -1
    //     sourcePatch = patchID
    // processor boundary:
    //     sourceProc = proc (on owner side)
    //     sourceNewNbrProc = distribution of coupled cell
    //     sourceFace = face (on owner side)
    //     sourcePatch = -1
    // ?cyclic:
    // ?    sourceProc = proc
    // ?    sourceNewNbrProc = distribution of coupled cell
    // ?    sourceFace = face (on owner side)
    // ?    sourcePatch = patchID
    // processor-cyclic boundary:
    //     sourceProc = proc (on owner side)
    //     sourceNewNbrProc = distribution of coupled cell
    //     sourceFace = face (on owner side)
    //     sourcePatch = patchID

    labelList sourcePatch;
    labelList sourceFace;
    labelList sourceProc;
    labelList sourceNewNbrProc;
    labelList sourcePointMaster;
    getCouplingData
    (
        distribution,
        sourceFace,
        sourceProc,
        sourcePatch,
        sourceNewNbrProc,
        sourcePointMaster
    );


    // Remove meshPhi. Since this would otherwise disappear anyway
    // during topo changes and we have to guarantee that all the fields
    // can be sent.
    mesh_.clearOut();
    mesh_.resetMotion();

    // Get data to send. Make sure is synchronised
    const wordList volScalars(mesh_.names(volScalarField::typeName));
    checkEqualWordList("volScalarFields", volScalars);
    const wordList volVectors(mesh_.names(volVectorField::typeName));
    checkEqualWordList("volVectorFields", volVectors);
    const wordList volSphereTensors
    (
        mesh_.names(volSphericalTensorField::typeName)
    );
    checkEqualWordList("volSphericalTensorFields", volSphereTensors);
    const wordList volSymmTensors(mesh_.names(volSymmTensorField::typeName));
    checkEqualWordList("volSymmTensorFields", volSymmTensors);
    const wordList volTensors(mesh_.names(volTensorField::typeName));
    checkEqualWordList("volTensorField", volTensors);

    const wordList surfScalars(mesh_.names(surfaceScalarField::typeName));
    checkEqualWordList("surfaceScalarFields", surfScalars);
    const wordList surfVectors(mesh_.names(surfaceVectorField::typeName));
    checkEqualWordList("surfaceVectorFields", surfVectors);
    const wordList surfSphereTensors
    (
        mesh_.names(surfaceSphericalTensorField::typeName)
    );
    checkEqualWordList("surfaceSphericalTensorFields", surfSphereTensors);
    const wordList surfSymmTensors
    (
        mesh_.names(surfaceSymmTensorField::typeName)
    );
    checkEqualWordList("surfaceSymmTensorFields", surfSymmTensors);
    const wordList surfTensors(mesh_.names(surfaceTensorField::typeName));
    checkEqualWordList("surfaceTensorFields", surfTensors);

    typedef volScalarField::Internal dimScalType;
    const wordList dimScalars(mesh_.names(dimScalType::typeName));
    checkEqualWordList("volScalarField::Internal", dimScalars);

    typedef volVectorField::Internal dimVecType;
    const wordList dimVectors(mesh_.names(dimVecType::typeName));
    checkEqualWordList("volVectorField::Internal", dimVectors);

    typedef volSphericalTensorField::Internal dimSphereType;
    const wordList dimSphereTensors(mesh_.names(dimSphereType::typeName));
    checkEqualWordList
    (
        "volSphericalTensorField::Internal",
        dimSphereTensors
    );

    typedef volSymmTensorField::Internal dimSymmTensorType;
    const wordList dimSymmTensors(mesh_.names(dimSymmTensorType::typeName));
    checkEqualWordList
    (
        "volSymmTensorField::Internal",
        dimSymmTensors
    );

    typedef volTensorField::Internal dimTensorType;
    const wordList dimTensors(mesh_.names(dimTensorType::typeName));
    checkEqualWordList("volTensorField::Internal", dimTensors);



    // Find patch to temporarily put exposed and processor faces into.
    label oldInternalPatchi = findNonEmptyPatch();



    // Delete processor patches, starting from the back. Move all faces into
    // oldInternalPatchi.
    labelList repatchFaceMap;
    {
        autoPtr<mapPolyMesh> repatchMap = deleteProcPatches(oldInternalPatchi);

        // Store face map (only face ordering that changed)
        repatchFaceMap = repatchMap().faceMap();

        // Reorder all boundary face data (sourceProc, sourceFace etc.)
        labelList bFaceMap
        (
            SubList<label>
            (
                repatchMap().reverseFaceMap(),
                mesh_.nFaces() - mesh_.nInternalFaces(),
                mesh_.nInternalFaces()
            )
          - mesh_.nInternalFaces()
        );

        inplaceReorder(bFaceMap, sourceFace);
        inplaceReorder(bFaceMap, sourceProc);
        inplaceReorder(bFaceMap, sourcePatch);
        inplaceReorder(bFaceMap, sourceNewNbrProc);
    }



    // Print a bit.
    if (debug)
    {
        Pout<< nl << "MESH WITH PROC PATCHES DELETED:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<volScalarField>(mesh_);
        printFieldInfo<volVectorField>(mesh_);
        printFieldInfo<volSphericalTensorField>(mesh_);
        printFieldInfo<volSymmTensorField>(mesh_);
        printFieldInfo<volTensorField>(mesh_);
        printFieldInfo<surfaceScalarField>(mesh_);
        printFieldInfo<surfaceVectorField>(mesh_);
        printFieldInfo<surfaceSphericalTensorField>(mesh_);
        printFieldInfo<surfaceSymmTensorField>(mesh_);
        printFieldInfo<surfaceTensorField>(mesh_);
        Pout<< nl << endl;
    }



    // Maps from subsetted mesh (that is sent) back to original maps
    labelListList subCellMap(Pstream::nProcs());
    labelListList subFaceMap(Pstream::nProcs());
    labelListList subPointMap(Pstream::nProcs());
    labelListList subPatchMap(Pstream::nProcs());
    // Maps from subsetted mesh to reconstructed mesh
    labelListList constructCellMap(Pstream::nProcs());
    labelListList constructFaceMap(Pstream::nProcs());
    labelListList constructPointMap(Pstream::nProcs());
    labelListList constructPatchMap(Pstream::nProcs());




    // Find out schedule
    // ~~~~~~~~~~~~~~~~~

    labelListList nSendCells(Pstream::nProcs());
    nSendCells[Pstream::myProcNo()] = countCells(distribution);
    Pstream::gatherList(nSendCells);
    Pstream::scatterList(nSendCells);


    // Allocate buffers
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);


    // What to send to neighbouring domains
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    bool oldParRun = UPstream::parRun();
    UPstream::parRun() = false;

    forAll(nSendCells[Pstream::myProcNo()], recvProc)
    {
        if
        (
            recvProc != Pstream::myProcNo()
         && nSendCells[Pstream::myProcNo()][recvProc] > 0
        )
        {
            // Send to recvProc

            if (debug)
            {
                Pout<< nl
                    << "SUBSETTING FOR DOMAIN " << recvProc
                    << " cells to send:"
                    << nSendCells[Pstream::myProcNo()][recvProc]
                    << nl << endl;
            }

            // Pstream for sending mesh and fields
            // OPstream str(Pstream::commsTypes::blocking, recvProc);
            UOPstream str(recvProc, pBufs);

            // Mesh subsetting engine
            fvMeshSubset subsetter(mesh_);

            // Subset the cells of the current domain.
            subsetter.setLargeCellSubset
            (
                distribution,
                recvProc,
                oldInternalPatchi,  // oldInternalFaces patch
                false               // no parallel sync
            );

            subCellMap[recvProc] = subsetter.cellMap();
            subFaceMap[recvProc] = subsetter.faceFlipMap();
            inplaceRenumberWithFlip
            (
                repatchFaceMap,
                false,      // oldToNew has flip
                true,       // subFaceMap has flip
                subFaceMap[recvProc]
            );
            subPointMap[recvProc] = subsetter.pointMap();
            subPatchMap[recvProc] = subsetter.patchMap();


            // Subset the boundary fields (owner/neighbour/processor)
            labelList procSourceFace;
            labelList procSourceProc;
            labelList procSourcePatch;
            labelList procSourceNewNbrProc;
            labelList procSourcePointMaster;

            subsetCouplingData
            (
                subsetter.subMesh(),
                subsetter.pointMap(),       // from subMesh to mesh
                subsetter.faceMap(),        //      ,,      ,,
                subsetter.cellMap(),        //      ,,      ,,

                distribution,               // old mesh distribution
                mesh_.faceOwner(),          // old owner
                mesh_.faceNeighbour(),
                mesh_.nInternalFaces(),

                sourceFace,
                sourceProc,
                sourcePatch,
                sourceNewNbrProc,
                sourcePointMaster,

                procSourceFace,
                procSourceProc,
                procSourcePatch,
                procSourceNewNbrProc,
                procSourcePointMaster
            );


            // Send to neighbour
            sendMesh
            (
                recvProc,
                subsetter.subMesh(),

                pointZoneNames,
                faceZoneNames,
                cellZoneNames,

                procSourceFace,
                procSourceProc,
                procSourcePatch,
                procSourceNewNbrProc,
                procSourcePointMaster,

                str
            );

            // volFields
            sendFields<volScalarField>(recvProc, volScalars, subsetter, str);
            sendFields<volVectorField>(recvProc, volVectors, subsetter, str);
            sendFields<volSphericalTensorField>
            (
                recvProc,
                volSphereTensors,
                subsetter,
                str
            );
            sendFields<volSymmTensorField>
            (
                recvProc,
                volSymmTensors,
                subsetter,
                str
            );
            sendFields<volTensorField>(recvProc, volTensors, subsetter, str);

            // surfaceFields
            sendFields<surfaceScalarField>
            (
                recvProc,
                surfScalars,
                subsetter,
                str
            );
            sendFields<surfaceVectorField>
            (
                recvProc,
                surfVectors,
                subsetter,
                str
            );
            sendFields<surfaceSphericalTensorField>
            (
                recvProc,
                surfSphereTensors,
                subsetter,
                str
            );
            sendFields<surfaceSymmTensorField>
            (
                recvProc,
                surfSymmTensors,
                subsetter,
                str
            );
            sendFields<surfaceTensorField>
            (
                recvProc,
                surfTensors,
                subsetter,
                str
            );

            // dimensionedFields
            sendFields<volScalarField::Internal>
            (
                recvProc,
                dimScalars,
                subsetter,
                str
            );
            sendFields<volVectorField::Internal>
            (
                recvProc,
                dimVectors,
                subsetter,
                str
            );
            sendFields<volSphericalTensorField::Internal>
            (
                recvProc,
                dimSphereTensors,
                subsetter,
                str
            );
            sendFields<volSymmTensorField::Internal>
            (
                recvProc,
                dimSymmTensors,
                subsetter,
                str
            );
            sendFields<volTensorField::Internal>
            (
                recvProc,
                dimTensors,
                subsetter,
                str
            );
        }
    }


    UPstream::parRun() = oldParRun;


    // Start sending&receiving from buffers
    pBufs.finishedSends();


    // Subset the part that stays
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // Save old mesh maps before changing mesh
        const labelList oldFaceOwner(mesh_.faceOwner());
        const labelList oldFaceNeighbour(mesh_.faceNeighbour());
        const label oldInternalFaces = mesh_.nInternalFaces();

        // Remove cells.
        autoPtr<mapPolyMesh> subMap
        (
            doRemoveCells
            (
                select(false, distribution, Pstream::myProcNo()),
                oldInternalPatchi
            )
        );

        // Addressing from subsetted mesh
        subCellMap[Pstream::myProcNo()] = subMap().cellMap();
        subFaceMap[Pstream::myProcNo()] = renumber
        (
            repatchFaceMap,
            subMap().faceMap()
        );
        // Insert the sign bit from face flipping
        labelList& faceMap = subFaceMap[Pstream::myProcNo()];
        forAll(faceMap, faceI)
        {
            faceMap[faceI] += 1;
        }
        const labelHashSet& flip = subMap().flipFaceFlux();
        forAllConstIter(labelHashSet, flip, iter)
        {
            label faceI = iter.key();
            faceMap[faceI] = -faceMap[faceI];
        }
        subPointMap[Pstream::myProcNo()] = subMap().pointMap();
        subPatchMap[Pstream::myProcNo()] = identity(patches.size());

        // Initialize all addressing into current mesh
        constructCellMap[Pstream::myProcNo()] = identity(mesh_.nCells());
        constructFaceMap[Pstream::myProcNo()] = identity(mesh_.nFaces()) + 1;
        constructPointMap[Pstream::myProcNo()] = identity(mesh_.nPoints());
        constructPatchMap[Pstream::myProcNo()] = identity(patches.size());

        // Subset the mesh data: neighbourCell/neighbourProc
        // fields
        labelList domainSourceFace;
        labelList domainSourceProc;
        labelList domainSourcePatch;
        labelList domainSourceNewNbrProc;
        labelList domainSourcePointMaster;

        subsetCouplingData
        (
            mesh_,                          // new mesh
            subMap().pointMap(),            // from new to original mesh
            subMap().faceMap(),             // from new to original mesh
            subMap().cellMap(),

            distribution,                   // distribution before subsetting
            oldFaceOwner,                   // owner before subsetting
            oldFaceNeighbour,               // neighbour        ,,
            oldInternalFaces,               // nInternalFaces   ,,

            sourceFace,
            sourceProc,
            sourcePatch,
            sourceNewNbrProc,
            sourcePointMaster,

            domainSourceFace,
            domainSourceProc,
            domainSourcePatch,
            domainSourceNewNbrProc,
            domainSourcePointMaster
        );

        sourceFace.transfer(domainSourceFace);
        sourceProc.transfer(domainSourceProc);
        sourcePatch.transfer(domainSourcePatch);
        sourceNewNbrProc.transfer(domainSourceNewNbrProc);
        sourcePointMaster.transfer(domainSourcePointMaster);
    }


    // Print a bit.
    if (debug)
    {
        Pout<< nl << "STARTING MESH:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<volScalarField>(mesh_);
        printFieldInfo<volVectorField>(mesh_);
        printFieldInfo<volSphericalTensorField>(mesh_);
        printFieldInfo<volSymmTensorField>(mesh_);
        printFieldInfo<volTensorField>(mesh_);
        printFieldInfo<surfaceScalarField>(mesh_);
        printFieldInfo<surfaceVectorField>(mesh_);
        printFieldInfo<surfaceSphericalTensorField>(mesh_);
        printFieldInfo<surfaceSymmTensorField>(mesh_);
        printFieldInfo<surfaceTensorField>(mesh_);
        Pout<< nl << endl;
    }



    // Receive and add what was sent
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    oldParRun = UPstream::parRun();
    UPstream::parRun() = false;

    forAll(nSendCells, sendProc)
    {
        // Did processor sendProc send anything to me?
        if
        (
            sendProc != Pstream::myProcNo()
         && nSendCells[sendProc][Pstream::myProcNo()] > 0
        )
        {
            if (debug)
            {
                Pout<< nl
                    << "RECEIVING FROM DOMAIN " << sendProc
                    << " cells to receive:"
                    << nSendCells[sendProc][Pstream::myProcNo()]
                    << nl << endl;
            }


            // Pstream for receiving mesh and fields
            UIPstream str(sendProc, pBufs);


            // Receive from sendProc
            labelList domainSourceFace;
            labelList domainSourceProc;
            labelList domainSourcePatch;
            labelList domainSourceNewNbrProc;
            labelList domainSourcePointMaster;

            autoPtr<fvMesh> domainMeshPtr;

            PtrList<volScalarField> vsf;
            PtrList<volVectorField> vvf;
            PtrList<volSphericalTensorField> vsptf;
            PtrList<volSymmTensorField> vsytf;
            PtrList<volTensorField> vtf;

            PtrList<surfaceScalarField> ssf;
            PtrList<surfaceVectorField> svf;
            PtrList<surfaceSphericalTensorField> ssptf;
            PtrList<surfaceSymmTensorField> ssytf;
            PtrList<surfaceTensorField> stf;

            PtrList<volScalarField::Internal> dsf;
            PtrList<volVectorField::Internal> dvf;
            PtrList<volSphericalTensorField::Internal> dstf;
            PtrList<volSymmTensorField::Internal> dsytf;
            PtrList<volTensorField::Internal> dtf;


            // Opposite of sendMesh
            {
                domainMeshPtr = receiveMesh
                (
                    sendProc,
                    pointZoneNames,
                    faceZoneNames,
                    cellZoneNames,

                    const_cast<Time&>(mesh_.time()),
                    domainSourceFace,
                    domainSourceProc,
                    domainSourcePatch,
                    domainSourceNewNbrProc,
                    domainSourcePointMaster,
                    str
                );
                fvMesh& domainMesh = domainMeshPtr();
                // Force construction of various on mesh.
                //(void)domainMesh.globalData();


                // Receive fields. Read as single dictionary because
                // of problems reading consecutive fields from single stream.
                dictionary fieldDicts(str);

                // Vol fields
                receiveFields<volScalarField>
                (
                    sendProc,
                    volScalars,
                    domainMesh,
                    vsf,
                    fieldDicts.subDict(volScalarField::typeName)
                );
                receiveFields<volVectorField>
                (
                    sendProc,
                    volVectors,
                    domainMesh,
                    vvf,
                    fieldDicts.subDict(volVectorField::typeName)
                );
                receiveFields<volSphericalTensorField>
                (
                    sendProc,
                    volSphereTensors,
                    domainMesh,
                    vsptf,
                    fieldDicts.subDict(volSphericalTensorField::typeName)
                );
                receiveFields<volSymmTensorField>
                (
                    sendProc,
                    volSymmTensors,
                    domainMesh,
                    vsytf,
                    fieldDicts.subDict(volSymmTensorField::typeName)
                );
                receiveFields<volTensorField>
                (
                    sendProc,
                    volTensors,
                    domainMesh,
                    vtf,
                    fieldDicts.subDict(volTensorField::typeName)
                );

                // Surface fields
                receiveFields<surfaceScalarField>
                (
                    sendProc,
                    surfScalars,
                    domainMesh,
                    ssf,
                    fieldDicts.subDict(surfaceScalarField::typeName)
                );
                receiveFields<surfaceVectorField>
                (
                    sendProc,
                    surfVectors,
                    domainMesh,
                    svf,
                    fieldDicts.subDict(surfaceVectorField::typeName)
                );
                receiveFields<surfaceSphericalTensorField>
                (
                    sendProc,
                    surfSphereTensors,
                    domainMesh,
                    ssptf,
                    fieldDicts.subDict(surfaceSphericalTensorField::typeName)
                );
                receiveFields<surfaceSymmTensorField>
                (
                    sendProc,
                    surfSymmTensors,
                    domainMesh,
                    ssytf,
                    fieldDicts.subDict(surfaceSymmTensorField::typeName)
                );
                receiveFields<surfaceTensorField>
                (
                    sendProc,
                    surfTensors,
                    domainMesh,
                    stf,
                    fieldDicts.subDict(surfaceTensorField::typeName)
                );

                // Dimensioned fields
                receiveFields<volScalarField::Internal>
                (
                    sendProc,
                    dimScalars,
                    domainMesh,
                    dsf,
                    fieldDicts.subDict
                    (
                        volScalarField::Internal::typeName
                    )
                );
                receiveFields<volVectorField::Internal>
                (
                    sendProc,
                    dimVectors,
                    domainMesh,
                    dvf,
                    fieldDicts.subDict
                    (
                        volVectorField::Internal::typeName
                    )
                );
                receiveFields<volSphericalTensorField::Internal>
                (
                    sendProc,
                    dimSphereTensors,
                    domainMesh,
                    dstf,
                    fieldDicts.subDict
                    (
                        volSphericalTensorField::Internal::
                        typeName
                    )
                );
                receiveFields<volSymmTensorField::Internal>
                (
                    sendProc,
                    dimSymmTensors,
                    domainMesh,
                    dsytf,
                    fieldDicts.subDict
                    (
                        volSymmTensorField::Internal::typeName
                    )
                );
                receiveFields<volTensorField::Internal>
                (
                    sendProc,
                    dimTensors,
                    domainMesh,
                    dtf,
                    fieldDicts.subDict
                    (
                        volTensorField::Internal::typeName
                    )
                );
            }
            const fvMesh& domainMesh = domainMeshPtr();


            constructCellMap[sendProc] = identity(domainMesh.nCells());
            constructFaceMap[sendProc] = identity(domainMesh.nFaces()) + 1;
            constructPointMap[sendProc] = identity(domainMesh.nPoints());
            constructPatchMap[sendProc] =
                identity(domainMesh.boundaryMesh().size());


            // Print a bit.
            if (debug)
            {
                Pout<< nl << "RECEIVED MESH FROM:" << sendProc << endl;
                printMeshInfo(domainMesh);
                printFieldInfo<volScalarField>(domainMesh);
                printFieldInfo<volVectorField>(domainMesh);
                printFieldInfo<volSphericalTensorField>(domainMesh);
                printFieldInfo<volSymmTensorField>(domainMesh);
                printFieldInfo<volTensorField>(domainMesh);
                printFieldInfo<surfaceScalarField>(domainMesh);
                printFieldInfo<surfaceVectorField>(domainMesh);
                printFieldInfo<surfaceSphericalTensorField>(domainMesh);
                printFieldInfo<surfaceSymmTensorField>(domainMesh);
                printFieldInfo<surfaceTensorField>(domainMesh);
            }


            // Now this mesh we received (from sendProc) needs to be merged
            // with the current mesh. On the current mesh we have for all
            // boundaryfaces the original face and processor. See if we can
            // match these up to the received domainSourceFace and
            // domainSourceProc.
            labelList masterCoupledFaces;
            labelList slaveCoupledFaces;
            findCouples
            (
                mesh_,

                sourceFace,
                sourceProc,
                sourcePatch,

                sendProc,
                domainMesh,
                domainSourceFace,
                domainSourceProc,
                domainSourcePatch,

                masterCoupledFaces,
                slaveCoupledFaces
            );

            // Generate additional coupling info (points, edges) from
            // faces-that-match
            faceCoupleInfo couples
            (
                mesh_,
                masterCoupledFaces,
                domainMesh,
                slaveCoupledFaces,
                mergeTol_,              // merge tolerance
                true,                   // faces align
                true,                   // couples are ordered already
                false
            );


            // Add domainMesh to mesh
            // ~~~~~~~~~~~~~~~~~~~~~~

            autoPtr<mapAddedPolyMesh> map = fvMeshAdder::add
            (
                mesh_,
                domainMesh,
                couples,
                false           // no parallel comms
            );

            // Update mesh data: sourceFace,sourceProc for added
            // mesh.

            sourceFace = mapBoundaryData
            (
                mesh_,
                map(),
                sourceFace,
                domainMesh.nInternalFaces(),
                domainSourceFace
            );
            sourceProc = mapBoundaryData
            (
                mesh_,
                map(),
                sourceProc,
                domainMesh.nInternalFaces(),
                domainSourceProc
            );
            sourcePatch = mapBoundaryData
            (
                mesh_,
                map(),
                sourcePatch,
                domainMesh.nInternalFaces(),
                domainSourcePatch
            );
            sourceNewNbrProc = mapBoundaryData
            (
                mesh_,
                map(),
                sourceNewNbrProc,
                domainMesh.nInternalFaces(),
                domainSourceNewNbrProc
            );
            // Update pointMaster data
            sourcePointMaster = mapPointData
            (
                mesh_,
                map(),
                sourcePointMaster,
                domainSourcePointMaster
            );


            // Update all addressing so xxProcAddressing points to correct
            // item in masterMesh.
            const labelList& oldCellMap = map().oldCellMap();
            const labelList& oldFaceMap = map().oldFaceMap();
            const labelList& oldPointMap = map().oldPointMap();
            const labelList& oldPatchMap = map().oldPatchMap();

            // Note: old mesh faces never flipped!
            forAll(constructPatchMap, proci)
            {
                if (proci != sendProc && constructPatchMap[proci].size())
                {
                    // Processor already in mesh (either myProcNo or received)
                    inplaceRenumber(oldCellMap, constructCellMap[proci]);
                    inplaceRenumberWithFlip
                    (
                        oldFaceMap,
                        false,
                        true,
                        constructFaceMap[proci]
                    );
                    inplaceRenumber(oldPointMap, constructPointMap[proci]);
                    inplaceRenumber(oldPatchMap, constructPatchMap[proci]);
                }
            }


            labelHashSet flippedAddedFaces;
            {
                // Find out if any faces of domain mesh were flipped (boundary
                // faces becoming internal)
                label nBnd = domainMesh.nFaces()-domainMesh.nInternalFaces();
                flippedAddedFaces.resize(nBnd/4);

                for
                (
                    label domainFaceI = domainMesh.nInternalFaces();
                    domainFaceI < domainMesh.nFaces();
                    domainFaceI++
                )
                {
                    label newFaceI = map().addedFaceMap()[domainFaceI];
                    label newCellI = mesh_.faceOwner()[newFaceI];

                    label domainCellI = domainMesh.faceOwner()[domainFaceI];

                    if (newCellI != map().addedCellMap()[domainCellI])
                    {
                        flippedAddedFaces.insert(domainFaceI);
                    }
                }
            }


            // Added processor
            inplaceRenumber(map().addedCellMap(), constructCellMap[sendProc]);
            // Add flip
            forAllConstIter(labelHashSet, flippedAddedFaces, iter)
            {
                label domainFaceI = iter.key();
                label& val = constructFaceMap[sendProc][domainFaceI];
                val = -val;
            }
            inplaceRenumberWithFlip
            (
                map().addedFaceMap(),
                false,
                true,           // constructFaceMap has flip sign
                constructFaceMap[sendProc]
            );
            inplaceRenumber(map().addedPointMap(), constructPointMap[sendProc]);
            inplaceRenumber(map().addedPatchMap(), constructPatchMap[sendProc]);

            if (debug)
            {
                Pout<< nl << "MERGED MESH FROM:" << sendProc << endl;
                printMeshInfo(mesh_);
                printFieldInfo<volScalarField>(mesh_);
                printFieldInfo<volVectorField>(mesh_);
                printFieldInfo<volSphericalTensorField>(mesh_);
                printFieldInfo<volSymmTensorField>(mesh_);
                printFieldInfo<volTensorField>(mesh_);
                printFieldInfo<surfaceScalarField>(mesh_);
                printFieldInfo<surfaceVectorField>(mesh_);
                printFieldInfo<surfaceSphericalTensorField>(mesh_);
                printFieldInfo<surfaceSymmTensorField>(mesh_);
                printFieldInfo<surfaceTensorField>(mesh_);
                Pout<< nl << endl;
            }
        }
    }

    UPstream::parRun() = oldParRun;

    // Print a bit.
    if (debug)
    {
        Pout<< nl << "REDISTRIBUTED MESH:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<volScalarField>(mesh_);
        printFieldInfo<volVectorField>(mesh_);
        printFieldInfo<volSphericalTensorField>(mesh_);
        printFieldInfo<volSymmTensorField>(mesh_);
        printFieldInfo<volTensorField>(mesh_);
        printFieldInfo<surfaceScalarField>(mesh_);
        printFieldInfo<surfaceVectorField>(mesh_);
        printFieldInfo<surfaceSphericalTensorField>(mesh_);
        printFieldInfo<surfaceSymmTensorField>(mesh_);
        printFieldInfo<surfaceTensorField>(mesh_);
        Pout<< nl << endl;
    }


    // See if any originally shared points need to be merged. Note: does
    // parallel comms. After this points and edges should again be consistent.
    mergeSharedPoints(sourcePointMaster, constructPointMap);


    // Add processorPatches
    // ~~~~~~~~~~~~~~~~~~~~

    // Per neighbour processor, per originating patch, the patchID
    // For faces resulting from internal faces or normal processor patches
    // the originating patch is -1. For cyclics this is the cyclic patchID.
    List<Map<label>> procPatchID;

    // Add processor and processorCyclic patches.
    addProcPatches(sourceNewNbrProc, sourcePatch, procPatchID);

    // Put faces into correct patch. Note that we now have proper
    // processorPolyPatches again so repatching will take care of coupled face
    // ordering.

    // Get boundary faces to be repatched. Is -1 or new patchID
    labelList newPatchID
    (
        getBoundaryPatch
        (
            sourceNewNbrProc,
            sourcePatch,
            procPatchID
        )
    );

    // Change patches. Since this might change ordering of coupled faces
    // we also need to adapt our constructMaps.
    repatch(newPatchID, constructFaceMap);

    // Bit of hack: processorFvPatchField does not get reset since created
    // from nothing so explicitly reset.
    initPatchFields<volScalarField, processorFvPatchField<scalar>>
    (
        Zero
    );
    initPatchFields<volVectorField, processorFvPatchField<vector>>
    (
        Zero
    );
    initPatchFields
    <
        volSphericalTensorField,
        processorFvPatchField<sphericalTensor>
    >
    (
        Zero
    );
    initPatchFields<volSymmTensorField, processorFvPatchField<symmTensor>>
    (
        Zero
    );
    initPatchFields<volTensorField, processorFvPatchField<tensor>>
    (
        Zero
    );


    mesh_.setInstance(mesh_.time().timeName());


    // Print a bit
    if (debug)
    {
        Pout<< nl << "FINAL MESH:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<volScalarField>(mesh_);
        printFieldInfo<volVectorField>(mesh_);
        printFieldInfo<volSphericalTensorField>(mesh_);
        printFieldInfo<volSymmTensorField>(mesh_);
        printFieldInfo<volTensorField>(mesh_);
        printFieldInfo<surfaceScalarField>(mesh_);
        printFieldInfo<surfaceVectorField>(mesh_);
        printFieldInfo<surfaceSphericalTensorField>(mesh_);
        printFieldInfo<surfaceSymmTensorField>(mesh_);
        printFieldInfo<surfaceTensorField>(mesh_);
        Pout<< nl << endl;
    }

    // Collect all maps and return
    return autoPtr<mapDistributePolyMesh>
    (
        new mapDistributePolyMesh
        (
            mesh_,

            nOldPoints,
            nOldFaces,
            nOldCells,
            oldPatchStarts.xfer(),
            oldPatchNMeshPoints.xfer(),

            subPointMap.xfer(),
            subFaceMap.xfer(),
            subCellMap.xfer(),
            subPatchMap.xfer(),

            constructPointMap.xfer(),
            constructFaceMap.xfer(),
            constructCellMap.xfer(),
            constructPatchMap.xfer(),

            true,           // subFaceMap has flip
            true            // constructFaceMap has flip
        )
    );
}


// ************************************************************************* //
