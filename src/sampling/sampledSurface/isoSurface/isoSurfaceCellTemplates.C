/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "isoSurfaceCell.H"
#include "polyMesh.H"
#include "tetMatcher.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Type Foam::isoSurfaceCell::generatePoint
(
    const DynamicList<Type>& snappedPoints,

    const scalar s0,
    const Type& p0,
    const label p0Index,

    const scalar s1,
    const Type& p1,
    const label p1Index
) const
{
    scalar d = s1-s0;

    if (mag(d) > VSMALL)
    {
        scalar s = (iso_-s0)/d;

        if (s >= 0.5 && s <= 1 && p1Index != -1)
        {
            return snappedPoints[p1Index];
        }
        else if (s >= 0.0 && s <= 0.5 && p0Index != -1)
        {
            return snappedPoints[p0Index];
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
void Foam::isoSurfaceCell::generateTriPoints
(
    const DynamicList<Type>& snapped,

    const scalar isoVal0,
    const scalar s0,
    const Type& p0,
    const label p0Index,

    const scalar isoVal1,
    const scalar s1,
    const Type& p1,
    const label p1Index,

    const scalar isoVal2,
    const scalar s2,
    const Type& p2,
    const label p2Index,

    const scalar isoVal3,
    const scalar s3,
    const Type& p3,
    const label p3Index,

    DynamicList<Type>& pts
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

        case 0x0E:
        case 0x01:
        {
            // 0 is common point. Orient such that normal points in positive
            // gradient direction
            if (isoVal0 >= isoVal1)
            {
                pts.append(generatePoint(snapped,s0,p0,p0Index,s1,p1,p1Index));
                pts.append(generatePoint(snapped,s0,p0,p0Index,s2,p2,p2Index));
                pts.append(generatePoint(snapped,s0,p0,p0Index,s3,p3,p3Index));
            }
            else
            {
                pts.append(generatePoint(snapped,s0,p0,p0Index,s2,p2,p2Index));
                pts.append(generatePoint(snapped,s0,p0,p0Index,s1,p1,p1Index));
                pts.append(generatePoint(snapped,s0,p0,p0Index,s3,p3,p3Index));
            }
        }
        break;

        case 0x0D:
        case 0x02:
        {
            // 1 is common point
            if (isoVal1 >= isoVal0)
            {
                pts.append(generatePoint(snapped,s1,p1,p1Index,s0,p0,p0Index));
                pts.append(generatePoint(snapped,s1,p1,p1Index,s3,p3,p3Index));
                pts.append(generatePoint(snapped,s1,p1,p1Index,s2,p2,p2Index));
            }
            else
            {
                pts.append(generatePoint(snapped,s1,p1,p1Index,s3,p3,p3Index));
                pts.append(generatePoint(snapped,s1,p1,p1Index,s0,p0,p0Index));
                pts.append(generatePoint(snapped,s1,p1,p1Index,s2,p2,p2Index));
            }
        }
        break;

        case 0x0C:
        case 0x03:
        {
            Type s02 = generatePoint(snapped,s0,p0,p0Index,s2,p2,p2Index);
            Type s13 = generatePoint(snapped,s1,p1,p1Index,s3,p3,p3Index);

            if (isoVal0 >= isoVal3)
            {
                pts.append(generatePoint(snapped,s0,p0,p0Index,s3,p3,p3Index));
                pts.append(s02);
                pts.append(s13);
                pts.append(s13);
                pts.append(generatePoint(snapped,s1,p1,p1Index,s2,p2,p2Index));
                pts.append(s02);
            }
            else
            {
                pts.append(s02);
                pts.append(generatePoint(snapped,s0,p0,p0Index,s3,p3,p3Index));
                pts.append(s13);
                pts.append(generatePoint(snapped,s1,p1,p1Index,s2,p2,p2Index));
                pts.append(s13);
                pts.append(s02);
            }
        }
        break;

        case 0x0B:
        case 0x04:
        {
            // 2 is common point
            if (isoVal2 >= isoVal0)
            {
                pts.append(generatePoint(snapped,s2,p2,p2Index,s0,p0,p0Index));
                pts.append(generatePoint(snapped,s2,p2,p2Index,s1,p1,p1Index));
                pts.append(generatePoint(snapped,s2,p2,p2Index,s3,p3,p3Index));
            }
            else
            {
                pts.append(generatePoint(snapped,s2,p2,p2Index,s1,p1,p1Index));
                pts.append(generatePoint(snapped,s2,p2,p2Index,s0,p0,p0Index));
                pts.append(generatePoint(snapped,s2,p2,p2Index,s3,p3,p3Index));
            }
        }
        break;

        case 0x0A:
        case 0x05:
        {
            Type s01 = generatePoint(snapped,s0,p0,p0Index,s1,p1,p1Index);
            Type s23 = generatePoint(snapped,s2,p2,p2Index,s3,p3,p3Index);

            if (isoVal3 >= isoVal0)
            {
                pts.append(s01);
                pts.append(s23);
                pts.append(generatePoint(snapped,s0,p0,p0Index,s3,p3,p3Index));
                pts.append(s01);
                pts.append(generatePoint(snapped,s1,p1,p1Index,s2,p2,p2Index));
                pts.append(s23);
            }
            else
            {
                pts.append(s23);
                pts.append(s01);
                pts.append(generatePoint(snapped,s0,p0,p0Index,s3,p3,p3Index));
                pts.append(generatePoint(snapped,s1,p1,p1Index,s2,p2,p2Index));
                pts.append(s01);
                pts.append(s23);
            }
        }
        break;

        case 0x09:
        case 0x06:
        {
            Type s01 = generatePoint(snapped,s0,p0,p0Index,s1,p1,p1Index);
            Type s23 = generatePoint(snapped,s2,p2,p2Index,s3,p3,p3Index);

            if (isoVal3 >= isoVal1)
            {
                pts.append(s01);
                pts.append(generatePoint(snapped,s1,p1,p1Index,s3,p3,p3Index));
                pts.append(s23);
                pts.append(s01);
                pts.append(generatePoint(snapped,s0,p0,p0Index,s2,p2,p2Index));
                pts.append(s23);
            }
            else
            {
                pts.append(generatePoint(snapped,s1,p1,p1Index,s3,p3,p3Index));
                pts.append(s01);
                pts.append(s23);
                pts.append(generatePoint(snapped,s0,p0,p0Index,s2,p2,p2Index));
                pts.append(s01);
                pts.append(s23);
            }
        }
        break;

        case 0x07:
        case 0x08:
        {
            // 3 is common point
            if (isoVal3 >= isoVal0)
            {
                pts.append(generatePoint(snapped,s3,p3,p3Index,s0,p0,p0Index));
                pts.append(generatePoint(snapped,s3,p3,p3Index,s2,p2,p2Index));
                pts.append(generatePoint(snapped,s3,p3,p3Index,s1,p1,p1Index));
            }
            else
            {
                pts.append(generatePoint(snapped,s3,p3,p3Index,s2,p2,p2Index));
                pts.append(generatePoint(snapped,s3,p3,p3Index,s0,p0,p0Index));
                pts.append(generatePoint(snapped,s3,p3,p3Index,s1,p1,p1Index));
            }
        }
        break;
    }
}


template<class Type>
void Foam::isoSurfaceCell::generateTriPoints
(
    const scalarField& cVals,
    const scalarField& pVals,

    const Field<Type>& cCoords,
    const Field<Type>& pCoords,

    const DynamicList<Type>& snappedPoints,
    const labelList& snappedCc,
    const labelList& snappedPoint,

    DynamicList<Type>& triPoints,
    DynamicList<label>& triMeshCells
) const
{
    tetMatcher tet;

    forAll(mesh_.cells(), cellI)
    {
        if (cellCutType_[cellI] != NOTCUT)
        {
            label oldNPoints = triPoints.size();

            const cell& cFaces = mesh_.cells()[cellI];

            if (tet.isA(mesh_, cellI))
            {
                // For tets don't do cell-centre decomposition, just use the
                // tet points and values

                const face& f0 = mesh_.faces()[cFaces[0]];

                // Get the other point
                const face& f1 = mesh_.faces()[cFaces[1]];
                label oppositeI = -1;
                forAll(f1, fp)
                {
                    oppositeI = f1[fp];

                    if (findIndex(f0, oppositeI) == -1)
                    {
                        break;
                    }
                }

                // Start off from positive volume tet to make sure we
                // generate outwards pointing tets
                if (mesh_.faceOwner()[cFaces[0]] == cellI)
                {
                    generateTriPoints
                    (
                        snappedPoints,

                        pVals_[f0[1]],
                        pVals[f0[1]],
                        pCoords[f0[1]],
                        snappedPoint[f0[1]],

                        pVals_[f0[0]],
                        pVals[f0[0]],
                        pCoords[f0[0]],
                        snappedPoint[f0[0]],

                        pVals_[f0[2]],
                        pVals[f0[2]],
                        pCoords[f0[2]],
                        snappedPoint[f0[2]],

                        pVals_[oppositeI],
                        pVals[oppositeI],
                        pCoords[oppositeI],
                        snappedPoint[oppositeI],

                        triPoints
                    );
                }
                else
                {
                    generateTriPoints
                    (
                        snappedPoints,

                        pVals_[f0[0]],
                        pVals[f0[0]],
                        pCoords[f0[0]],
                        snappedPoint[f0[0]],

                        pVals_[f0[1]],
                        pVals[f0[1]],
                        pCoords[f0[1]],
                        snappedPoint[f0[1]],

                        pVals_[f0[2]],
                        pVals[f0[2]],
                        pCoords[f0[2]],
                        snappedPoint[f0[2]],

                        pVals_[oppositeI],
                        pVals[oppositeI],
                        pCoords[oppositeI],
                        snappedPoint[oppositeI],

                        triPoints
                    );
                }
            }
            else
            {
                const cell& cFaces = mesh_.cells()[cellI];

                forAll(cFaces, cFaceI)
                {
                    label faceI = cFaces[cFaceI];
                    const face& f = mesh_.faces()[faceI];

                    const label fp0 = mesh_.tetBasePtIs()[faceI];

                    label fp = f.fcIndex(fp0);
                    for (label i = 2; i < f.size(); i++)
                    {
                        label nextFp = f.fcIndex(fp);
                        triFace tri(f[fp0], f[fp], f[nextFp]);

                        // Start off from positive volume tet to make sure we
                        // generate outwards pointing tets
                        if (mesh_.faceOwner()[faceI] == cellI)
                        {
                            generateTriPoints
                            (
                                snappedPoints,

                                pVals_[tri[1]],
                                pVals[tri[1]],
                                pCoords[tri[1]],
                                snappedPoint[tri[1]],

                                pVals_[tri[0]],
                                pVals[tri[0]],
                                pCoords[tri[0]],
                                snappedPoint[tri[0]],

                                pVals_[tri[2]],
                                pVals[tri[2]],
                                pCoords[tri[2]],
                                snappedPoint[tri[2]],

                                cVals_[cellI],
                                cVals[cellI],
                                cCoords[cellI],
                                snappedCc[cellI],

                                triPoints
                            );
                        }
                        else
                        {
                            generateTriPoints
                            (
                                snappedPoints,

                                pVals_[tri[0]],
                                pVals[tri[0]],
                                pCoords[tri[0]],
                                snappedPoint[tri[0]],

                                pVals_[tri[1]],
                                pVals[tri[1]],
                                pCoords[tri[1]],
                                snappedPoint[tri[1]],

                                pVals_[tri[2]],
                                pVals[tri[2]],
                                pCoords[tri[2]],
                                snappedPoint[tri[2]],

                                cVals_[cellI],
                                cVals[cellI],
                                cCoords[cellI],
                                snappedCc[cellI],

                                triPoints
                            );
                        }

                        fp = nextFp;
                    }
                }
            }


            // Every three triPoints is a cell
            label nCells = (triPoints.size()-oldNPoints)/3;
            for (label i = 0; i < nCells; i++)
            {
                triMeshCells.append(cellI);
            }
        }
    }

    triPoints.shrink();
    triMeshCells.shrink();
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::isoSurfaceCell::interpolate
(
    const scalarField& cVals,
    const scalarField& pVals,
    const Field<Type>& cCoords,
    const Field<Type>& pCoords
) const
{
    DynamicList<Type> triPoints(nCutCells_);
    DynamicList<label> triMeshCells(nCutCells_);

    // Dummy snap data
    DynamicList<Type> snappedPoints;
    labelList snappedCc(mesh_.nCells(), -1);
    labelList snappedPoint(mesh_.nPoints(), -1);


    generateTriPoints
    (
        cVals,
        pVals,

        cCoords,
        pCoords,

        snappedPoints,
        snappedCc,
        snappedPoint,

        triPoints,
        triMeshCells
    );


    // One value per point
    tmp<Field<Type> > tvalues(new Field<Type>(points().size()));
    Field<Type>& values = tvalues();

    forAll(triPoints, i)
    {
        label mergedPointI = triPointMergeMap_[i];

        if (mergedPointI >= 0)
        {
            values[mergedPointI] = triPoints[i];
        }
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::isoSurfaceCell::interpolate
(
    const Field<Type>& cCoords,
    const Field<Type>& pCoords
) const
{
    DynamicList<Type> triPoints(nCutCells_);
    DynamicList<label> triMeshCells(nCutCells_);

    // Dummy snap data
    DynamicList<Type> snappedPoints;
    labelList snappedCc(mesh_.nCells(), -1);
    labelList snappedPoint(mesh_.nPoints(), -1);


    generateTriPoints
    (
        cVals_,
        pVals_,

        cCoords,
        pCoords,

        snappedPoints,
        snappedCc,
        snappedPoint,

        triPoints,
        triMeshCells
    );


    // One value per point
    tmp<Field<Type> > tvalues(new Field<Type>(points().size()));
    Field<Type>& values = tvalues();

    forAll(triPoints, i)
    {
        label mergedPointI = triPointMergeMap_[i];

        if (mergedPointI >= 0)
        {
            values[mergedPointI] = triPoints[i];
        }
    }

    return tvalues;
}


// ************************************************************************* //
