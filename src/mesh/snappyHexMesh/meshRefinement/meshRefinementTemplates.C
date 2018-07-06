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
#include "fvMesh.H"
#include "globalIndex.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Add a T entry
template<class T> void Foam::meshRefinement::updateList
(
    const labelList& newToOld,
    const T& nullValue,
    List<T>& elems
)
{
    List<T> newElems(newToOld.size(), nullValue);

    forAll(newElems, i)
    {
        label oldI = newToOld[i];

        if (oldI >= 0)
        {
            newElems[i] = elems[oldI];
        }
    }

    elems.transfer(newElems);
}


template<class T>
T Foam::meshRefinement::gAverage
(
    const PackedBoolList& isMasterElem,
    const UList<T>& values
)
{
    if (values.size() != isMasterElem.size())
    {
        FatalErrorInFunction
            << "Number of elements in list " << values.size()
            << " does not correspond to number of elements in isMasterElem "
            << isMasterElem.size()
            << exit(FatalError);
    }

    T sum = Zero;
    label n = 0;

    forAll(values, i)
    {
        if (isMasterElem[i])
        {
            sum += values[i];
            n++;
        }
    }

    reduce(sum, sumOp<T>());
    reduce(n, sumOp<label>());

    if (n > 0)
    {
        return sum/n;
    }
    else
    {
        return pTraits<T>::max;
    }
}


// Compare two lists over all boundary faces
template<class T>
void Foam::meshRefinement::testSyncBoundaryFaceList
(
    const scalar tol,
    const string& msg,
    const UList<T>& faceData,
    const UList<T>& syncedFaceData
) const
{
    label nBFaces = mesh_.nFaces() - mesh_.nInternalFaces();

    if (faceData.size() != nBFaces || syncedFaceData.size() != nBFaces)
    {
        FatalErrorInFunction
            << "Boundary faces:" << nBFaces
            << " faceData:" << faceData.size()
            << " syncedFaceData:" << syncedFaceData.size()
            << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        label bFacei = pp.start() - mesh_.nInternalFaces();

        forAll(pp, i)
        {
            const T& data = faceData[bFacei];
            const T& syncData = syncedFaceData[bFacei];

            if (mag(data - syncData) > tol)
            {
                label facei = pp.start()+i;

                FatalErrorInFunction
                    << msg
                    << "patchFace:" << i
                    << " face:" << facei
                    << " fc:" << mesh_.faceCentres()[facei]
                    << " patch:" << pp.name()
                    << " faceData:" << data
                    << " syncedFaceData:" << syncData
                    << " diff:" << mag(data - syncData)
                    << abort(FatalError);
            }

            bFacei++;
        }
    }
}


// Print list sorted by coordinates. Used for comparing non-parallel v.s.
// parallel operation
template<class T>
void Foam::meshRefinement::collectAndPrint
(
    const UList<point>& points,
    const UList<T>& data
)
{
    globalIndex globalPoints(points.size());

    pointField allPoints;
    globalPoints.gather
    (
        Pstream::worldComm,
        identity(Pstream::nProcs()),
        points,
        allPoints,
        UPstream::msgType(),
        Pstream::commsTypes::blocking
    );

    List<T> allData;
    globalPoints.gather
    (
        Pstream::worldComm,
        identity(Pstream::nProcs()),
        data,
        allData,
        UPstream::msgType(),
        Pstream::commsTypes::blocking
    );


    scalarField magAllPoints(mag(allPoints-point(-0.317, 0.117, 0.501)));

    labelList visitOrder;
    sortedOrder(magAllPoints, visitOrder);
    forAll(visitOrder, i)
    {
        label allPointi = visitOrder[i];
        Info<< allPoints[allPointi] << " : " << allData[allPointi]
            << endl;
    }
}


template<class GeoField>
void Foam::meshRefinement::addPatchFields
(
    fvMesh& mesh,
    const word& patchFieldType
)
{
    HashTable<GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        GeoField& fld = *iter();
        typename GeoField::Boundary& fldBf =
            fld.boundaryFieldRef();

        label sz = fldBf.size();
        fldBf.setSize(sz+1);
        fldBf.set
        (
            sz,
            GeoField::Patch::New
            (
                patchFieldType,
                mesh.boundary()[sz],
                fld()
            )
        );
    }
}


template<class GeoField>
void Foam::meshRefinement::reorderPatchFields
(
    fvMesh& mesh,
    const labelList& oldToNew
)
{
    HashTable<GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        iter()->boundaryFieldRef().reorder(oldToNew);
    }
}


template<class Enum>
int Foam::meshRefinement::readFlags
(
    const Enum& namedEnum,
    const wordList& words
)
{
    int flags = 0;

    forAll(words, i)
    {
        int index = namedEnum[words[i]];
        int val = 1<<index;
        flags |= val;
    }
    return flags;
}


template<class Type>
void Foam::meshRefinement::weightedSum
(
    const polyMesh& mesh,
    const PackedBoolList& isMasterEdge,
    const labelList& meshPoints,
    const edgeList& edges,
    const scalarField& edgeWeights,
    const Field<Type>& pointData,
    Field<Type>& sum
)
{
    if
    (
        edges.size() != isMasterEdge.size()
     || edges.size() != edgeWeights.size()
     || meshPoints.size() != pointData.size()
    )
    {
        FatalErrorInFunction
            << "Inconsistent sizes for edge or point data:"
            << " isMasterEdge:" << isMasterEdge.size()
            << " edgeWeights:" << edgeWeights.size()
            << " edges:" << edges.size()
            << " pointData:" << pointData.size()
            << " meshPoints:" << meshPoints.size()
            << abort(FatalError);
    }

    sum.setSize(meshPoints.size());
    sum = Zero;

    forAll(edges, edgeI)
    {
        if (isMasterEdge[edgeI])
        {
            const edge& e = edges[edgeI];

            scalar eWeight = edgeWeights[edgeI];

            label v0 = e[0];
            label v1 = e[1];

            sum[v0] += eWeight*pointData[v1];
            sum[v1] += eWeight*pointData[v0];
        }
    }

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        sum,
        plusEqOp<Type>(),
        Type(Zero)     // null value
    );
}


// ************************************************************************* //
