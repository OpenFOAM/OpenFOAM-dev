/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
#include "fvMeshAdder.H"
#include "processorCyclicFvPatchField.H"
#include "polyTopoChange.H"
#include "removeCells.H"
#include "polyModifyFace.H"
#include "polyDistributionMap.H"
#include "syncTools.H"
#include "CompactListList.H"
#include "fvMeshTools.H"
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
    const labelList& referNbrPatchID_;

public:

    lessProcPatches
    (
        const labelList& nbrProc,
        const labelList& referPatchID,
        const labelList& referNbrPatchID
    )
    :
        nbrProc_(nbrProc),
        referPatchID_(referPatchID),
        referNbrPatchID_(referNbrPatchID)
    {}

    bool operator()(const label a, const label b)
    {
        // Lower processor ID-s go first
        if (nbrProc_[a] < nbrProc_[b])
        {
            return true;
        }
        else if (nbrProc_[a] > nbrProc_[b])
        {
            return false;
        }

        // Non-cyclics go next
        else if (referPatchID_[a] == -1)
        {
            return true;
        }
        else if (referPatchID_[b] == -1)
        {
            return false;
        }

        // Cyclics should be ordered by refer patch ID if this is the owner
        // (lower processor ID), and by the neighbour refer patch ID if this is
        // the neighbour
        else if (Pstream::myProcNo() < nbrProc_[a])
        {
            return referPatchID_[a] < referPatchID_[b];
        }
        else
        {
            return referNbrPatchID_[a] < referNbrPatchID_[b];
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


Foam::wordList Foam::fvMeshDistribute::fieldNames
(
    const word& typeName,
    label& nFields
) const
{
    wordList fieldNames(mesh_.names(typeName));

    if (fieldNames.size())
    {
        HashSet<word> fieldSet(fieldNames);
        fieldSet -= fvMesh::geometryFields;
        fieldNames = fieldSet.toc();
        nFields += checkEqualWordList(typeName, fieldNames);
    }

    return fieldNames;
}


Foam::label Foam::fvMeshDistribute::checkEqualWordList
(
    const word& typeName,
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
                << "When checking for equal numbers of " << typeName
                << " :" << nl
                << "processor0 has:" << allNames[0] << nl
                << "processor" << proci << " has:" << allNames[proci] << nl
                << typeName << " need to be synchronised on all processors."
                << exit(FatalError);
        }
    }

    return lst.size();
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


Foam::label Foam::fvMeshDistribute::findInternalPatch() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label internalPatchi = -1;

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (isA<internalPolyPatch>(pp))
        {
            internalPatchi = patchi;
            break;
        }
    }

    if (internalPatchi == -1)
    {
        FatalErrorInFunction
            << "Cannot find a internal patch in " << patches.names() << nl
            << "    of types " << patches.types() << nl
            << "    An internal patch must be provided for the exposed "
               "internal faces." << exit(FatalError);
    }

    if (debug)
    {
        Pout<< "findInternalPatch : using patch " << internalPatchi
            << " name:" << patches[internalPatchi].name()
            << " type:" << patches[internalPatchi].type()
            << " for the exposed internal faces." << endl;
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

    return internalPatchi;
}

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
            << "Cannot find a non-empty patch in " << patches.names() << nl
            << "    of types " << patches.types() << nl
            << "    An non-empty patch must be provided for the exposed "
               "internal faces." << exit(FatalError);
    }

    if (debug)
    {
        Pout<< "findNonEmptyPatch : using patch " << nonEmptyPatchi
            << " name:" << patches[nonEmptyPatchi].name()
            << " type:" << patches[nonEmptyPatchi].type()
            << " for the exposed non-empty faces." << endl;
    }

    return nonEmptyPatchi;
}


Foam::autoPtr<Foam::polyTopoChangeMap> Foam::fvMeshDistribute::deleteProcPatches
(
    const label destinationPatch
)
{
    // New patchID per boundary faces to be repatched. Is -1 (no change)
    // or new patchID
    labelList newPatchID(mesh_.nFaces() - mesh_.nInternalFaces(), -1);

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
        }
    }

    // Note: order of boundary faces been kept the same since the
    // destinationPatch is at the end and we have visited the patches in
    // incremental order.
    labelListList dummyFaceMaps;
    autoPtr<polyTopoChangeMap> map = repatch(newPatchID, dummyFaceMaps);


    // Delete (now empty) processor patches.
    {
        labelList oldToNew(identityMap(mesh_.boundaryMesh().size()));
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


Foam::autoPtr<Foam::polyTopoChangeMap> Foam::fvMeshDistribute::repatch
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
    // is currently not supported by topoChange.

    // Store boundary fields (we only do this for surfaceFields)
    PtrList<FieldField<fvsPatchField, scalar>> sFields;
    saveBoundaryFields<scalar, surfaceMesh>(sFields);
    PtrList<FieldField<fvsPatchField, vector>> vFields;
    saveBoundaryFields<vector, surfaceMesh>(vFields);
    PtrList<FieldField<fvsPatchField, sphericalTensor>> sptFields;
    saveBoundaryFields<sphericalTensor, surfaceMesh>(sptFields);
    PtrList<FieldField<fvsPatchField, symmTensor>> sytFields;
    saveBoundaryFields<symmTensor, surfaceMesh>(sytFields);
    PtrList<FieldField<fvsPatchField, tensor>> tFields;
    saveBoundaryFields<tensor, surfaceMesh>(tFields);

    // Change the mesh (no inflation). Note: parallel comms allowed.
    //
    // NOTE: there is one very particular problem with this ordering.
    // We first create the processor patches and use these to merge out
    // shared points (see mergeSharedPoints below). So temporarily points
    // and edges do not match!

    autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.mapFields(map);

    // Map patch fields using stored boundary fields. Note: assumes order
    // of fields has not changed in object registry!
    mapBoundaryFields<scalar, surfaceMesh>(map, sFields);
    mapBoundaryFields<vector, surfaceMesh>(map, vFields);
    mapBoundaryFields<sphericalTensor, surfaceMesh>(map, sptFields);
    mapBoundaryFields<symmTensor, surfaceMesh>(map, sytFields);
    mapBoundaryFields<tensor, surfaceMesh>(map, tFields);

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


Foam::autoPtr<Foam::polyTopoChangeMap> Foam::fvMeshDistribute::mergeSharedPoints
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
        return autoPtr<polyTopoChangeMap>(nullptr);
    }

    // Create the mesh change engine to merge the points
    polyTopoChange meshMod(mesh_);
    {
        // Remove all non-master points.
        forAll(mesh_.points(), pointi)
        {
            Map<label>::const_iterator iter = pointToMaster.find(pointi);

            if (iter != pointToMaster.end())
            {
                if (iter() != pointi)
                {
                    meshMod.removePoint(pointi, iter());
                }
            }
        }

        // Modify faces for points. Note: could use pointFaces here but want to
        // avoid addressing calculation.
        const faceList& faces = mesh_.faces();

        forAll(faces, facei)
        {
            const face& f = faces[facei];

            bool hasMerged = false;

            forAll(f, fp)
            {
                label pointi = f[fp];

                Map<label>::const_iterator iter = pointToMaster.find(pointi);

                if (iter != pointToMaster.end())
                {
                    if (iter() != pointi)
                    {
                        hasMerged = true;
                        break;
                    }
                }
            }

            if (hasMerged)
            {
                face newF(f);

                forAll(f, fp)
                {
                    label pointi = f[fp];

                    Map<label>::const_iterator iter =
                        pointToMaster.find(pointi);

                    if (iter != pointToMaster.end())
                    {
                        newF[fp] = iter();
                    }
                }

                label patchID = mesh_.boundaryMesh().whichPatch(facei);
                label nei = (patchID == -1 ? mesh_.faceNeighbour()[facei] : -1);
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
                        newF,                       // modified face
                        facei,                      // label of face
                        mesh_.faceOwner()[facei],   // owner
                        nei,                        // neighbour
                        false,                      // face flip
                        patchID,                    // patch for face
                        false,                      // remove from zone
                        zoneID,                     // zone for face
                        zoneFlip                    // face flip in zone
                    )
                );
            }
        }
    }

    // Change the mesh (no inflation). Note: parallel comms allowed.
    autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.mapFields(map);

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
    labelList& sourceNbrPatch,
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
    sourceNbrPatch.setSize(nBnd);
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


            label patchi = -1, nbrPatchi = -1;
            if (isA<processorCyclicPolyPatch>(pp))
            {
                patchi =
                    refCast<const processorCyclicPolyPatch>(pp)
                   .referPatchID();
                nbrPatchi =
                    refCast<const cyclicPolyPatch>(patches[patchi])
                   .nbrPatchID();

            }

            forAll(pp, i)
            {
                label bndI = offset + i;
                sourcePatch[bndI] = patchi;
                sourceNbrPatch[bndI] = nbrPatchi;
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
                    sourceNbrPatch[bndI] = cpp.nbrPatchID();
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
                    sourceNbrPatch[bndI] = cpp.nbrPatchID();
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
                sourceNbrPatch[bndI] = -1;
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
        const distributionMap& slavesMap = gmd.globalCoPointSlavesMap();
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
    const labelList& sourceNbrPatch,
    const labelList& sourceNewNbrProc,
    const labelList& sourcePointMaster,

    labelList& subFace,
    labelList& subProc,
    labelList& subPatch,
    labelList& subNbrPatch,
    labelList& subNewNbrProc,
    labelList& subPointMaster
)
{
    subFace.setSize(mesh.nFaces() - mesh.nInternalFaces());
    subProc.setSize(mesh.nFaces() - mesh.nInternalFaces());
    subPatch.setSize(mesh.nFaces() - mesh.nInternalFaces());
    subNbrPatch.setSize(mesh.nFaces() - mesh.nInternalFaces());
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
            subNbrPatch[newBFacei] = sourceNbrPatch[oldBFacei];
            subNewNbrProc[newBFacei] = sourceNewNbrProc[oldBFacei];
        }
    }


    subPointMaster = UIndirectList<label>(sourcePointMaster, pointMap);
}


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


Foam::autoPtr<Foam::polyTopoChangeMap> Foam::fvMeshDistribute::doRemoveCells
(
    const labelList& cellsToRemove,
    const label oldInternalPatchi
)
{
    // Mesh change engine
    polyTopoChange meshMod(mesh_);

    // Cell removal topo engine. Do NOT synchronise parallel since
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

    // Save surface fields. This is not done as GeometricField as these would
    // get mapped. Fields are flattened for convenience.
    PtrList<Field<scalar>> sFields;
    PtrList<Field<vector>> vFields;
    PtrList<Field<sphericalTensor>> sptFields;
    PtrList<Field<symmTensor>> sytFields;
    PtrList<Field<tensor>> tFields;
    initMapExposedFaces(sFields);
    initMapExposedFaces(vFields);
    initMapExposedFaces(sptFields);
    initMapExposedFaces(sytFields);
    initMapExposedFaces(tFields);

    // Change the mesh. No inflation. Note: no parallel comms allowed.
    autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh_, false, false);

    // Update fields
    mesh_.mapFields(map);

    // Any exposed faces in a surfaceField will not be mapped. Map the value
    // of these separately (until there is support in all PatchFields for
    // mapping from internal faces ...)
    mapExposedFaces(map(), sFields);
    mapExposedFaces(map(), vFields);
    mapExposedFaces(map(), sptFields);
    mapExposedFaces(map(), sytFields);
    mapExposedFaces(map(), tFields);

    return map;
}


void Foam::fvMeshDistribute::addProcPatches
(
    const labelList& nbrProc,         // Processor that neighbour is now on
    const labelList& referPatchID,    // Original patch ID (or -1)
    const labelList& referNbrPatchID, // Original neighbour patch ID (or -1)
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
    sortedOrder
    (
        nbrProc,
        indices,
        lessProcPatches(nbrProc, referPatchID, referNbrPatchID)
    );

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
                        pcPatch.name()
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
    const labelList& sourceNbrPatch,
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
        const meshPointZones& pointZones = mesh.pointZones();

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
        const meshFaceZones& faceZones = mesh.faceZones();

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
        const meshCellZones& cellZones = mesh.cellZones();

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
    //    const meshCellZones& cellZones = mesh.cellZones();
    //
    //    forAll(cellZones, zoneI)
    //    {
    //        UIndirectList<label>(cellZoneID, cellZones[zoneI]) = zoneI;
    //    }
    //}

    // Send
    toDomain
        << mesh.points()
        << CompactListList<label>(mesh.faces())
        << mesh.faceOwner()
        << mesh.faceNeighbour()
        << mesh.boundaryMesh()

        //*** Write the old-time volumes if present
        // << mesh.V0().field()
        // << mesh.V0().field()

        << zonePoints
        << zoneFaces
        << zoneFaceFlip
        << zoneCells

        << sourceFace
        << sourceProc
        << sourcePatch
        << sourceNbrPatch
        << sourceNewNbrProc
        << sourcePointMaster;


    if (debug)
    {
        Pout<< "Started sending mesh to domain " << domain
            << endl;
    }
}


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
    labelList& domainSourceNbrPatch,
    labelList& domainSourceNewNbrProc,
    labelList& domainSourcePointMaster,
    Istream& fromNbr
)
{
    pointField domainPoints(fromNbr);
    faceList domainFaces = CompactListList<label>(fromNbr).list<face>();
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
        >> domainSourceNbrPatch
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
                runTime.name(),
                runTime,
                IOobject::NO_READ
            ),
            move(domainPoints),
            move(domainFaces),
            move(domainAllOwner),
            move(domainAllNeighbour),
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

Foam::fvMeshDistribute::fvMeshDistribute(fvMesh& mesh)
:
    mesh_(mesh)
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


Foam::autoPtr<Foam::polyDistributionMap> Foam::fvMeshDistribute::distribute
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
        return autoPtr<polyDistributionMap>
        (
            new polyDistributionMap
            (
                mesh_,

                nOldPoints,
                nOldFaces,
                nOldCells,
                move(oldPatchStarts),
                move(oldPatchNMeshPoints),

                labelListList(1, identityMap(mesh_.nPoints())),
                labelListList(1, identityMap(mesh_.nFaces())),
                labelListList(1, identityMap(mesh_.nCells())),
                labelListList(1, identityMap(patches.size())),

                labelListList(1, identityMap(mesh_.nPoints())),
                labelListList(1, identityMap(mesh_.nFaces())),
                labelListList(1, identityMap(mesh_.nCells())),
                labelListList(1, identityMap(patches.size()))
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

    labelList sourceFace;
    labelList sourceProc;
    labelList sourcePatch;
    labelList sourceNbrPatch;
    labelList sourceNewNbrProc;
    labelList sourcePointMaster;
    getCouplingData
    (
        distribution,
        sourceFace,
        sourceProc,
        sourcePatch,
        sourceNbrPatch,
        sourceNewNbrProc,
        sourcePointMaster
    );


    // Remove old-time geometry to avoid the need to distribute it
    mesh_.resetMotion();

    label nFields = 0;

    // Get data to send. Make sure is synchronised
    const wordList volScalars
    (
        fieldNames(volScalarField::typeName, nFields)
    );
    const wordList volVectors
    (
        fieldNames(volVectorField::typeName, nFields)
    );
    const wordList volSphereTensors
    (
        fieldNames(volSphericalTensorField::typeName, nFields)
    );
    const wordList volSymmTensors(fieldNames
    (
        volSymmTensorField::typeName, nFields)
    );
    const wordList volTensors
    (
        fieldNames(volTensorField::typeName, nFields)
    );

    const wordList surfScalars
    (
        fieldNames(surfaceScalarField::typeName, nFields)
    );
    const wordList surfVectors
    (
        fieldNames(surfaceVectorField::typeName, nFields)
    );
    const wordList surfSphereTensors
    (
        fieldNames(surfaceSphericalTensorField::typeName, nFields)
    );
    const wordList surfSymmTensors
    (
        fieldNames(surfaceSymmTensorField::typeName, nFields)
    );
    const wordList surfTensors
    (
        fieldNames(surfaceTensorField::typeName, nFields)
    );

    const wordList pointScalars
    (
        fieldNames(pointScalarField::typeName, nFields)
    );
    const wordList pointVectors
    (
        fieldNames(pointVectorField::typeName, nFields)
    );
    const wordList pointSphereTensors
    (
        fieldNames(pointSphericalTensorField::typeName, nFields)
    );
    const wordList pointSymmTensors
    (
        fieldNames(pointSymmTensorField::typeName, nFields)
    );
    const wordList pointTensors
    (
        fieldNames(pointTensorField::typeName, nFields)
    );

    const wordList dimScalars
    (
        fieldNames(volScalarField::Internal::typeName, nFields)
    );
    const wordList dimVectors
    (
        fieldNames(volVectorField::Internal::typeName, nFields)
    );
    const wordList dimSphereTensors
    (
        fieldNames(volSphericalTensorField::Internal::typeName, nFields)
    );
    const wordList dimSymmTensors
    (
        fieldNames(volSymmTensorField::Internal::typeName, nFields)
    );
    const wordList dimTensors
    (
        fieldNames(volTensorField::Internal::typeName, nFields)
    );

    // Find patch to temporarily put exposed internal and processor faces into.
    // If there are no fields patch 0 is used,
    // If there are fields the internal patch is used.
    const label oldInternalPatchi =
        nFields
      ? findInternalPatch()
      : findNonEmptyPatch();

    // Delete processor patches, starting from the back. Move all faces into
    // oldInternalPatchi.
    labelList repatchFaceMap;
    {
        autoPtr<polyTopoChangeMap> repatchMap =
            deleteProcPatches(oldInternalPatchi);

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
        inplaceReorder(bFaceMap, sourceNbrPatch);
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
        printFieldInfo<pointScalarField>(mesh_);
        printFieldInfo<pointVectorField>(mesh_);
        printFieldInfo<pointSphericalTensorField>(mesh_);
        printFieldInfo<pointSymmTensorField>(mesh_);
        printFieldInfo<pointTensorField>(mesh_);
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
            labelList procSourceNbrPatch;
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
                sourceNbrPatch,
                sourceNewNbrProc,
                sourcePointMaster,

                procSourceFace,
                procSourceProc,
                procSourcePatch,
                procSourceNbrPatch,
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
                procSourceNbrPatch,
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

            // pointFields
            sendFields<pointScalarField>
            (
                recvProc,
                pointScalars,
                subsetter,
                str
            );
            sendFields<pointVectorField>
            (
                recvProc,
                pointVectors,
                subsetter,
                str
            );
            sendFields<pointSphericalTensorField>
            (
                recvProc,
                pointSphereTensors,
                subsetter,
                str
            );
            sendFields<pointSymmTensorField>
            (
                recvProc,
                pointSymmTensors,
                subsetter,
                str
            );
            sendFields<pointTensorField>
            (
                recvProc,
                pointTensors,
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
        autoPtr<polyTopoChangeMap> subMap
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
        forAll(faceMap, facei)
        {
            faceMap[facei] += 1;
        }
        const labelHashSet& flip = subMap().flipFaceFlux();
        forAllConstIter(labelHashSet, flip, iter)
        {
            label facei = iter.key();
            faceMap[facei] = -faceMap[facei];
        }
        subPointMap[Pstream::myProcNo()] = subMap().pointMap();
        subPatchMap[Pstream::myProcNo()] = identityMap(patches.size());

        // Initialise all addressing into current mesh
        constructCellMap[Pstream::myProcNo()] = identityMap(mesh_.nCells());
        constructFaceMap[Pstream::myProcNo()] = identityMap(mesh_.nFaces()) + 1;
        constructPointMap[Pstream::myProcNo()] = identityMap(mesh_.nPoints());
        constructPatchMap[Pstream::myProcNo()] = identityMap(patches.size());

        // Subset the mesh data: neighbourCell/neighbourProc
        // fields
        labelList domainSourceFace;
        labelList domainSourceProc;
        labelList domainSourcePatch;
        labelList domainSourceNbrPatch;
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
            sourceNbrPatch,
            sourceNewNbrProc,
            sourcePointMaster,

            domainSourceFace,
            domainSourceProc,
            domainSourcePatch,
            domainSourceNbrPatch,
            domainSourceNewNbrProc,
            domainSourcePointMaster
        );

        sourceFace.transfer(domainSourceFace);
        sourceProc.transfer(domainSourceProc);
        sourcePatch.transfer(domainSourcePatch);
        sourceNbrPatch.transfer(domainSourceNbrPatch);
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
        printFieldInfo<pointScalarField>(mesh_);
        printFieldInfo<pointVectorField>(mesh_);
        printFieldInfo<pointSphericalTensorField>(mesh_);
        printFieldInfo<pointSymmTensorField>(mesh_);
        printFieldInfo<pointTensorField>(mesh_);
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
            labelList domainSourceNbrPatch;
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

            PtrList<pointScalarField> psf;
            PtrList<pointVectorField> pvf;
            PtrList<pointSphericalTensorField> psptf;
            PtrList<pointSymmTensorField> psytf;
            PtrList<pointTensorField> ptf;

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
                    domainSourceNbrPatch,
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

                // Point fields
                pointMesh& domainPointMesh =
                    const_cast<pointMesh&>(pointMesh::New(domainMesh));
                receiveFields<pointScalarField>
                (
                    sendProc,
                    pointScalars,
                    domainPointMesh,
                    psf,
                    fieldDicts.subDict(pointScalarField::typeName)
                );
                receiveFields<pointVectorField>
                (
                    sendProc,
                    pointVectors,
                    domainPointMesh,
                    pvf,
                    fieldDicts.subDict(pointVectorField::typeName)
                );
                receiveFields<pointSphericalTensorField>
                (
                    sendProc,
                    pointSphereTensors,
                    domainPointMesh,
                    psptf,
                    fieldDicts.subDict(pointSphericalTensorField::typeName)
                );
                receiveFields<pointSymmTensorField>
                (
                    sendProc,
                    pointSymmTensors,
                    domainPointMesh,
                    psytf,
                    fieldDicts.subDict(pointSymmTensorField::typeName)
                );
                receiveFields<pointTensorField>
                (
                    sendProc,
                    pointTensors,
                    domainPointMesh,
                    ptf,
                    fieldDicts.subDict(pointTensorField::typeName)
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


            constructCellMap[sendProc] = identityMap(domainMesh.nCells());
            constructFaceMap[sendProc] = identityMap(domainMesh.nFaces()) + 1;
            constructPointMap[sendProc] = identityMap(domainMesh.nPoints());
            constructPatchMap[sendProc] =
                identityMap(domainMesh.boundaryMesh().size());


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
                printFieldInfo<pointScalarField>(domainMesh);
                printFieldInfo<pointVectorField>(domainMesh);
                printFieldInfo<pointSphericalTensorField>(domainMesh);
                printFieldInfo<pointSymmTensorField>(domainMesh);
                printFieldInfo<pointTensorField>(domainMesh);
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
                slaveCoupledFaces
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
            sourceNbrPatch = mapBoundaryData
            (
                mesh_,
                map(),
                sourceNbrPatch,
                domainMesh.nInternalFaces(),
                domainSourceNbrPatch
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
                    label domainFacei = domainMesh.nInternalFaces();
                    domainFacei < domainMesh.nFaces();
                    domainFacei++
                )
                {
                    label newFacei = map().addedFaceMap()[domainFacei];
                    label newCellI = mesh_.faceOwner()[newFacei];

                    label domainCellI = domainMesh.faceOwner()[domainFacei];

                    if (newCellI != map().addedCellMap()[domainCellI])
                    {
                        flippedAddedFaces.insert(domainFacei);
                    }
                }
            }


            // Added processor
            inplaceRenumber(map().addedCellMap(), constructCellMap[sendProc]);
            // Add flip
            forAllConstIter(labelHashSet, flippedAddedFaces, iter)
            {
                label domainFacei = iter.key();
                label& val = constructFaceMap[sendProc][domainFacei];
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
                printFieldInfo<pointScalarField>(mesh_);
                printFieldInfo<pointVectorField>(mesh_);
                printFieldInfo<pointSphericalTensorField>(mesh_);
                printFieldInfo<pointSymmTensorField>(mesh_);
                printFieldInfo<pointTensorField>(mesh_);
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
        printFieldInfo<pointScalarField>(mesh_);
        printFieldInfo<pointVectorField>(mesh_);
        printFieldInfo<pointSphericalTensorField>(mesh_);
        printFieldInfo<pointSymmTensorField>(mesh_);
        printFieldInfo<pointTensorField>(mesh_);
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
    addProcPatches(sourceNewNbrProc, sourcePatch, sourceNbrPatch, procPatchID);

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

    // Correct coupled patch fields
    correctCoupledPatchFields<volScalarField>();
    correctCoupledPatchFields<volVectorField>();
    correctCoupledPatchFields<volSphericalTensorField>();
    correctCoupledPatchFields<volSymmTensorField>();
    correctCoupledPatchFields<volTensorField>();

    mesh_.setInstance(mesh_.time().name());

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
        printFieldInfo<pointScalarField>(mesh_);
        printFieldInfo<pointVectorField>(mesh_);
        printFieldInfo<pointSphericalTensorField>(mesh_);
        printFieldInfo<pointSymmTensorField>(mesh_);
        printFieldInfo<pointTensorField>(mesh_);
        Pout<< nl << endl;
    }

    // Collect all maps and return
    return autoPtr<polyDistributionMap>
    (
        new polyDistributionMap
        (
            mesh_,

            nOldPoints,
            nOldFaces,
            nOldCells,
            move(oldPatchStarts),
            move(oldPatchNMeshPoints),

            move(subPointMap),
            move(subFaceMap),
            move(subCellMap),
            move(subPatchMap),

            move(constructPointMap),
            move(constructFaceMap),
            move(constructCellMap),
            move(constructPatchMap),

            true,           // subFaceMap has flip
            true            // constructFaceMap has flip
        )
    );
}


// ************************************************************************* //
