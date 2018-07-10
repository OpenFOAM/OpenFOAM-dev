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

#include "isoSurface.H"
#include "polyMesh.H"
#include "syncTools.H"
#include "surfaceFields.H"
#include "OFstream.H"
#include "meshTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::SlicedGeometricField
<
    Type,
    Foam::fvPatchField,
    Foam::slicedFvPatchField,
    Foam::volMesh
>>
Foam::isoSurface::adaptPatchFields
(
    const GeometricField<Type, fvPatchField, volMesh>& fld
) const
{
    typedef SlicedGeometricField
    <
        Type,
        fvPatchField,
        slicedFvPatchField,
        volMesh
    > FieldType;

    tmp<FieldType> tsliceFld
    (
        new FieldType
        (
            IOobject
            (
                fld.name(),
                fld.instance(),
                fld.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fld,        // internal field
            true        // preserveCouples
        )
    );
    FieldType& sliceFld = tsliceFld.ref();

    const fvMesh& mesh = fld.mesh();

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    typename FieldType::Boundary& sliceFldBf =
        sliceFld.boundaryFieldRef();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if
        (
            isA<emptyPolyPatch>(pp)
         && pp.size() != sliceFldBf[patchi].size()
        )
        {
            // Clear old value. Cannot resize it since is a slice.
            sliceFldBf.set(patchi, nullptr);

            // Set new value we can change
            sliceFldBf.set
            (
                patchi,
                new calculatedFvPatchField<Type>
                (
                    mesh.boundary()[patchi],
                    sliceFld
                )
            );

            // Note: cannot use patchInternalField since uses emptyFvPatch::size
            // Do our own internalField instead.
            const labelUList& faceCells =
                mesh.boundary()[patchi].patch().faceCells();

            Field<Type>& pfld = sliceFldBf[patchi];
            pfld.setSize(faceCells.size());
            forAll(faceCells, i)
            {
                pfld[i] = sliceFld[faceCells[i]];
            }
        }
        else if (isA<cyclicPolyPatch>(pp))
        {
            // Already has interpolate as value
        }
        else if (isA<processorPolyPatch>(pp))
        {
            fvPatchField<Type>& pfld = const_cast<fvPatchField<Type>&>
            (
                sliceFldBf[patchi]
            );

            const scalarField& w = mesh.weights().boundaryField()[patchi];

            tmp<Field<Type>> f =
                w*pfld.patchInternalField()
              + (1.0-w)*pfld.patchNeighbourField();

            PackedBoolList isCollocated
            (
                collocatedFaces(refCast<const processorPolyPatch>(pp))
            );

            forAll(isCollocated, i)
            {
                if (!isCollocated[i])
                {
                    pfld[i] = f()[i];
                }
            }
        }
    }
    return tsliceFld;
}


template<class Type>
Type Foam::isoSurface::generatePoint
(
    const scalar s0,
    const Type& p0,
    const bool hasSnap0,
    const Type& snapP0,

    const scalar s1,
    const Type& p1,
    const bool hasSnap1,
    const Type& snapP1
) const
{
    scalar d = s1-s0;

    if (mag(d) > vSmall)
    {
        scalar s = (iso_-s0)/d;

        if (hasSnap1 && s >= 0.5 && s <= 1)
        {
            return snapP1;
        }
        else if (hasSnap0 && s >= 0.0 && s <= 0.5)
        {
            return snapP0;
        }
        else
        {
            return s*p1 + (1.0-s)*p0;
        }
    }
    else
    {
        scalar s = 0.4999;

        return s*p1 + (1.0-s)*p0;
    }
}


template<class Type>
void Foam::isoSurface::generateTriPoints
(
    const scalar s0,
    const Type& p0,
    const bool hasSnap0,
    const Type& snapP0,

    const scalar s1,
    const Type& p1,
    const bool hasSnap1,
    const Type& snapP1,

    const scalar s2,
    const Type& p2,
    const bool hasSnap2,
    const Type& snapP2,

    const scalar s3,
    const Type& p3,
    const bool hasSnap3,
    const Type& snapP3,

    DynamicList<Type>& points
) const
{
    int triIndex = 0;
    if (s0 < iso_)
    {
        triIndex |= 1;
    }
    if (s1 < iso_)
    {
        triIndex |= 2;
    }
    if (s2 < iso_)
    {
        triIndex |= 4;
    }
    if (s3 < iso_)
    {
        triIndex |= 8;
    }

    /* Form the vertices of the triangles for each case */
    switch (triIndex)
    {
        case 0x00:
        case 0x0F:
        break;

        case 0x01:
        case 0x0E:
            points.append
            (
                generatePoint(s0,p0,hasSnap0,snapP0,s1,p1,hasSnap1,snapP1)
            );
            points.append
            (
                generatePoint(s0,p0,hasSnap0,snapP0,s2,p2,hasSnap2,snapP2)
            );
            points.append
            (
                generatePoint(s0,p0,hasSnap0,snapP0,s3,p3,hasSnap3,snapP3)
            );
            if (triIndex == 0x0E)
            {
                // Flip normals
                label sz = points.size();
                Swap(points[sz-2], points[sz-1]);
            }
        break;

        case 0x02:
        case 0x0D:
            points.append
            (
                generatePoint(s1,p1,hasSnap1,snapP1,s0,p0,hasSnap0,snapP0)
            );
            points.append
            (
                generatePoint(s1,p1,hasSnap1,snapP1,s3,p3,hasSnap3,snapP3)
            );
            points.append
            (
                generatePoint(s1,p1,hasSnap1,snapP1,s2,p2,hasSnap2,snapP2)
            );
            if (triIndex == 0x0D)
            {
                // Flip normals
                label sz = points.size();
                Swap(points[sz-2], points[sz-1]);
            }
        break;

        case 0x03:
        case 0x0C:
            {
                Type p0p2 =
                    generatePoint(s0,p0,hasSnap0,snapP0,s2,p2,hasSnap2,snapP2);
                Type p1p3 =
                    generatePoint(s1,p1,hasSnap1,snapP1,s3,p3,hasSnap3,snapP3);

                points.append
                (
                    generatePoint(s0,p0,hasSnap0,snapP0,s3,p3,hasSnap3,snapP3)
                );
                points.append(p1p3);
                points.append(p0p2);

                points.append(p1p3);
                points.append
                (
                    generatePoint(s1,p1,hasSnap1,snapP1,s2,p2,hasSnap2,snapP2)
                );
                points.append(p0p2);
            }
            if (triIndex == 0x0C)
            {
                // Flip normals
                label sz = points.size();
                Swap(points[sz-5], points[sz-4]);
                Swap(points[sz-2], points[sz-1]);
            }
        break;

        case 0x04:
        case 0x0B:
            {
                points.append
                (
                    generatePoint(s2,p2,hasSnap2,snapP2,s0,p0,hasSnap0,snapP0)
                );
                points.append
                (
                    generatePoint(s2,p2,hasSnap2,snapP2,s1,p1,hasSnap1,snapP1)
                );
                points.append
                (
                    generatePoint(s2,p2,hasSnap2,snapP2,s3,p3,hasSnap3,snapP3)
                );
            }
            if (triIndex == 0x0B)
            {
                // Flip normals
                label sz = points.size();
                Swap(points[sz-2], points[sz-1]);
            }
        break;

        case 0x05:
        case 0x0A:
            {
                Type p0p1 =
                    generatePoint(s0,p0,hasSnap0,snapP0,s1,p1,hasSnap1,snapP1);
                Type p2p3 =
                    generatePoint(s2,p2,hasSnap2,snapP2,s3,p3,hasSnap3,snapP3);

                points.append(p0p1);
                points.append(p2p3);
                points.append
                (
                    generatePoint(s0,p0,hasSnap0,snapP0,s3,p3,hasSnap3,snapP3)
                );

                points.append(p0p1);
                points.append
                (
                    generatePoint(s1,p1,hasSnap1,snapP1,s2,p2,hasSnap2,snapP2)
                );
                points.append(p2p3);
            }
            if (triIndex == 0x0A)
            {
                // Flip normals
                label sz = points.size();
                Swap(points[sz-5], points[sz-4]);
                Swap(points[sz-2], points[sz-1]);
            }
        break;

        case 0x06:
        case 0x09:
            {
                Type p0p1 =
                    generatePoint(s0,p0,hasSnap0,snapP0,s1,p1,hasSnap1,snapP1);
                Type p2p3 =
                    generatePoint(s2,p2,hasSnap2,snapP2,s3,p3,hasSnap3,snapP3);

                points.append(p0p1);
                points.append
                (
                    generatePoint(s1,p1,hasSnap1,snapP1,s3,p3,hasSnap3,snapP3)
                );
                points.append(p2p3);

                points.append(p0p1);
                points.append(p2p3);
                points.append
                (
                    generatePoint(s0,p0,hasSnap0,snapP0,s2,p2,hasSnap2,snapP2)
                );
            }
            if (triIndex == 0x09)
            {
                // Flip normals
                label sz = points.size();
                Swap(points[sz-5], points[sz-4]);
                Swap(points[sz-2], points[sz-1]);
            }
        break;

        case 0x08:
        case 0x07:
            points.append
            (
                generatePoint(s3,p3,hasSnap3,snapP3,s0,p0,hasSnap0,snapP0)
            );
            points.append
            (
                generatePoint(s3,p3,hasSnap3,snapP3,s2,p2,hasSnap2,snapP2)
            );
            points.append
            (
                generatePoint(s3,p3,hasSnap3,snapP3,s1,p1,hasSnap1,snapP1)
            );
            if (triIndex == 0x07)
            {
                // Flip normals
                label sz = points.size();
                Swap(points[sz-2], points[sz-1]);
            }
        break;
    }
}


template<class Type>
Foam::label Foam::isoSurface::generateFaceTriPoints
(
    const volScalarField& cVals,
    const scalarField& pVals,

    const GeometricField<Type, fvPatchField, volMesh>& cCoords,
    const Field<Type>& pCoords,

    const DynamicList<Type>& snappedPoints,
    const labelList& snappedCc,
    const labelList& snappedPoint,
    const label facei,

    const scalar neiVal,
    const Type& neiPt,
    const bool hasNeiSnap,
    const Type& neiSnapPt,

    DynamicList<Type>& triPoints,
    DynamicList<label>& triMeshCells
) const
{
    label own = mesh_.faceOwner()[facei];

    label oldNPoints = triPoints.size();

    const face& f = mesh_.faces()[facei];

    forAll(f, fp)
    {
        label pointi = f[fp];
        label nextPointi = f[f.fcIndex(fp)];

        generateTriPoints
        (
            pVals[pointi],
            pCoords[pointi],
            snappedPoint[pointi] != -1,
            (
                snappedPoint[pointi] != -1
              ? snappedPoints[snappedPoint[pointi]]
              : Type(Zero)
            ),

            pVals[nextPointi],
            pCoords[nextPointi],
            snappedPoint[nextPointi] != -1,
            (
                snappedPoint[nextPointi] != -1
              ? snappedPoints[snappedPoint[nextPointi]]
              : Type(Zero)
            ),

            cVals[own],
            cCoords[own],
            snappedCc[own] != -1,
            (
                snappedCc[own] != -1
              ? snappedPoints[snappedCc[own]]
              : Type(Zero)
            ),

            neiVal,
            neiPt,
            hasNeiSnap,
            neiSnapPt,

            triPoints
        );
    }

    // Every three triPoints is a triangle
    label nTris = (triPoints.size()-oldNPoints)/3;
    for (label i = 0; i < nTris; i++)
    {
        triMeshCells.append(own);
    }

    return nTris;
}


template<class Type>
void Foam::isoSurface::generateTriPoints
(
    const volScalarField& cVals,
    const scalarField& pVals,

    const GeometricField<Type, fvPatchField, volMesh>& cCoords,
    const Field<Type>& pCoords,

    const DynamicList<Type>& snappedPoints,
    const labelList& snappedCc,
    const labelList& snappedPoint,

    DynamicList<Type>& triPoints,
    DynamicList<label>& triMeshCells
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    if
    (
        (cVals.size() != mesh_.nCells())
     || (pVals.size() != mesh_.nPoints())
     || (cCoords.size() != mesh_.nCells())
     || (pCoords.size() != mesh_.nPoints())
     || (snappedCc.size() != mesh_.nCells())
     || (snappedPoint.size() != mesh_.nPoints())
    )
    {
        FatalErrorInFunction
            << "Incorrect size." << endl
            << "mesh: nCells:" << mesh_.nCells()
            << " points:" << mesh_.nPoints() << endl
            << "cVals:" << cVals.size() << endl
            << "cCoords:" << cCoords.size() << endl
            << "snappedCc:" << snappedCc.size() << endl
            << "pVals:" << pVals.size() << endl
            << "pCoords:" << pCoords.size() << endl
            << "snappedPoint:" << snappedPoint.size() << endl
            << abort(FatalError);
    }


    // Generate triangle points

    triPoints.clear();
    triMeshCells.clear();

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (faceCutType_[facei] != NOTCUT)
        {
            generateFaceTriPoints
            (
                cVals,
                pVals,

                cCoords,
                pCoords,

                snappedPoints,
                snappedCc,
                snappedPoint,
                facei,

                cVals[nei[facei]],
                cCoords[nei[facei]],
                snappedCc[nei[facei]] != -1,
                (
                    snappedCc[nei[facei]] != -1
                  ? snappedPoints[snappedCc[nei[facei]]]
                  : Type(Zero)
                ),

                triPoints,
                triMeshCells
            );
        }
    }


    // Determine neighbouring snap status
    boolList neiSnapped(mesh_.nFaces()-mesh_.nInternalFaces(), false);
    List<Type> neiSnappedPoint(neiSnapped.size(), Type(Zero));
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            label facei = pp.start();
            forAll(pp, i)
            {
                label bFacei = facei-mesh_.nInternalFaces();
                label snappedIndex = snappedCc[own[facei]];

                if (snappedIndex != -1)
                {
                    neiSnapped[bFacei] = true;
                    neiSnappedPoint[bFacei] = snappedPoints[snappedIndex];
                }
                facei++;
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh_, neiSnapped);
    syncTools::swapBoundaryFaceList(mesh_, neiSnappedPoint);



    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& cpp =
                refCast<const processorPolyPatch>(pp);

            PackedBoolList isCollocated(collocatedFaces(cpp));

            forAll(isCollocated, i)
            {
                label facei = pp.start()+i;

                if (faceCutType_[facei] != NOTCUT)
                {
                    if (isCollocated[i])
                    {
                        generateFaceTriPoints
                        (
                            cVals,
                            pVals,

                            cCoords,
                            pCoords,

                            snappedPoints,
                            snappedCc,
                            snappedPoint,
                            facei,

                            cVals.boundaryField()[patchi][i],
                            cCoords.boundaryField()[patchi][i],
                            neiSnapped[facei-mesh_.nInternalFaces()],
                            neiSnappedPoint[facei-mesh_.nInternalFaces()],

                            triPoints,
                            triMeshCells
                        );
                    }
                    else
                    {
                        generateFaceTriPoints
                        (
                            cVals,
                            pVals,

                            cCoords,
                            pCoords,

                            snappedPoints,
                            snappedCc,
                            snappedPoint,
                            facei,

                            cVals.boundaryField()[patchi][i],
                            cCoords.boundaryField()[patchi][i],
                            false,
                            Type(Zero),

                            triPoints,
                            triMeshCells
                        );
                    }
                }
            }
        }
        else
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                if (faceCutType_[facei] != NOTCUT)
                {
                    generateFaceTriPoints
                    (
                        cVals,
                        pVals,

                        cCoords,
                        pCoords,

                        snappedPoints,
                        snappedCc,
                        snappedPoint,
                        facei,

                        cVals.boundaryField()[patchi][i],
                        cCoords.boundaryField()[patchi][i],
                        false,  // fc not snapped
                        Type(Zero),

                        triPoints,
                        triMeshCells
                    );
                }
                facei++;
            }
        }
    }

    triPoints.shrink();
    triMeshCells.shrink();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::isoSurface::interpolate
(
    const label nPoints,
    const labelList& triPointMergeMap,
    const DynamicList<Type>& unmergedValues
)
{
    // One value per point
    tmp<Field<Type>> tvalues(new Field<Type>(nPoints, Type(Zero)));
    Field<Type>& values = tvalues.ref();
    labelList nValues(values.size(), 0);

    forAll(unmergedValues, i)
    {
        label mergedPointi = triPointMergeMap[i];

        if (mergedPointi >= 0)
        {
            values[mergedPointi] += unmergedValues[i];
            nValues[mergedPointi]++;
        }
    }

    forAll(values, i)
    {
        if (nValues[i] > 0)
        {
            values[i] /= scalar(nValues[i]);
        }
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::isoSurface::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& cCoords,
    const Field<Type>& pCoords
) const
{
    // Recalculate boundary values
    tmp<SlicedGeometricField
    <
        Type,
        fvPatchField,
        slicedFvPatchField,
        volMesh
    >> c2(adaptPatchFields(cCoords));


    DynamicList<Type> triPoints(3*nCutCells_);
    DynamicList<label> triMeshCells(nCutCells_);

    // Dummy snap data
    DynamicList<Type> snappedPoints;
    labelList snappedCc(mesh_.nCells(), -1);
    labelList snappedPoint(mesh_.nPoints(), -1);

    generateTriPoints
    (
        cValsPtr_(),
        pVals_,

        c2(),
        pCoords,

        snappedPoints,
        snappedCc,
        snappedPoint,

        triPoints,
        triMeshCells
    );

    return interpolate
    (
        points().size(),
        triPointMergeMap_,
        triPoints
    );
}


// ************************************************************************* //
