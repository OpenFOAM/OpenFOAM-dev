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

#include "isoSurfaceCell.H"
#include "polyMesh.H"
#include "tetMatcher.H"
#include "isoSurface.H"

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

    if (mag(d) > vSmall)
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

    const scalar s0,
    const Type& p0,
    const label p0Index,

    const scalar s1,
    const Type& p1,
    const label p1Index,

    const scalar s2,
    const Type& p2,
    const label p2Index,

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

        case 0x01:
        case 0x0E:
        {
            pts.append
            (
                generatePoint(snapped,s0,p0,p0Index,s1,p1,p1Index)
            );
            pts.append
            (
                generatePoint(snapped,s0,p0,p0Index,s2,p2,p2Index)
            );
            pts.append
            (
                generatePoint(snapped,s0,p0,p0Index,s3,p3,p3Index)
            );
            if (triIndex == 0x0E)
            {
                // Flip normals
                label sz = pts.size();
                Swap(pts[sz-2], pts[sz-1]);
            }
        }
        break;

        case 0x02:
        case 0x0D:
        {
            pts.append
            (
                generatePoint(snapped,s1,p1,p1Index,s0,p0,p0Index)
            );
            pts.append
            (
                generatePoint(snapped,s1,p1,p1Index,s3,p3,p3Index)
            );
            pts.append
            (
                generatePoint(snapped,s1,p1,p1Index,s2,p2,p2Index)
            );

            if (triIndex == 0x0D)
            {
                // Flip normals
                label sz = pts.size();
                Swap(pts[sz-2], pts[sz-1]);
            }
        }
        break;

        case 0x03:
        case 0x0C:
        {
            Type p0p2 =
                generatePoint(snapped,s0,p0,p0Index,s2,p2,p2Index);
            Type p1p3 =
                generatePoint(snapped,s1,p1,p1Index,s3,p3,p3Index);

            pts.append
            (
                generatePoint(snapped,s0,p0,p0Index,s3,p3,p3Index)
            );
            pts.append(p1p3);
            pts.append(p0p2);

            pts.append(p1p3);
            pts.append
            (
                generatePoint(snapped,s1,p1,p1Index,s2,p2,p2Index)
            );
            pts.append(p0p2);

            if (triIndex == 0x0C)
            {
                // Flip normals
                label sz = pts.size();
                Swap(pts[sz-5], pts[sz-4]);
                Swap(pts[sz-2], pts[sz-1]);
            }
        }
        break;

        case 0x04:
        case 0x0B:
        {
            pts.append
            (
                generatePoint(snapped,s2,p2,p2Index,s0,p0,p0Index)
            );
            pts.append
            (
                generatePoint(snapped,s2,p2,p2Index,s1,p1,p1Index)
            );
            pts.append
            (
                generatePoint(snapped,s2,p2,p2Index,s3,p3,p3Index)
            );

            if (triIndex == 0x0B)
            {
                // Flip normals
                label sz = pts.size();
                Swap(pts[sz-2], pts[sz-1]);
            }
        }
        break;

        case 0x05:
        case 0x0A:
        {
            Type p0p1 =
                generatePoint(snapped,s0,p0,p0Index,s1,p1,p1Index);
            Type p2p3 =
                generatePoint(snapped,s2,p2,p2Index,s3,p3,p3Index);

            pts.append(p0p1);
            pts.append(p2p3);
            pts.append
            (
                generatePoint(snapped,s0,p0,p0Index,s3,p3,p3Index)
            );

            pts.append(p0p1);
            pts.append
            (
                generatePoint(snapped,s1,p1,p1Index,s2,p2,p2Index)
            );
            pts.append(p2p3);

            if (triIndex == 0x0A)
            {
                // Flip normals
                label sz = pts.size();
                Swap(pts[sz-5], pts[sz-4]);
                Swap(pts[sz-2], pts[sz-1]);
            }
        }
        break;

        case 0x06:
        case 0x09:
        {
            Type p0p1 =
                generatePoint(snapped,s0,p0,p0Index,s1,p1,p1Index);
            Type p2p3 =
                generatePoint(snapped,s2,p2,p2Index,s3,p3,p3Index);

            pts.append(p0p1);
            pts.append
            (
                generatePoint(snapped,s1,p1,p1Index,s3,p3,p3Index)
            );
            pts.append(p2p3);

            pts.append(p0p1);
            pts.append(p2p3);
            pts.append
            (
                generatePoint(snapped,s0,p0,p0Index,s2,p2,p2Index)
            );

            if (triIndex == 0x09)
            {
                // Flip normals
                label sz = pts.size();
                Swap(pts[sz-5], pts[sz-4]);
                Swap(pts[sz-2], pts[sz-1]);
            }
        }
        break;

        case 0x08:
        case 0x07:
        {
            pts.append
            (
                generatePoint(snapped,s3,p3,p3Index,s0,p0,p0Index)
            );
            pts.append
            (
                generatePoint(snapped,s3,p3,p3Index,s2,p2,p2Index)
            );
            pts.append
            (
                generatePoint(snapped,s3,p3,p3Index,s1,p1,p1Index)
            );
            if (triIndex == 0x07)
            {
                // Flip normals
                label sz = pts.size();
                Swap(pts[sz-2], pts[sz-1]);
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
    label countNotFoundTets = 0;

    forAll(mesh_.cells(), celli)
    {
        if (cellCutType_[celli] != NOTCUT)
        {
            label oldNPoints = triPoints.size();

            const cell& cFaces = mesh_.cells()[celli];

            if (tet.isA(mesh_, celli))
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
                if (mesh_.faceOwner()[cFaces[0]] == celli)
                {
                    generateTriPoints
                    (
                        snappedPoints,

                        pVals[f0[1]],
                        pCoords[f0[1]],
                        snappedPoint[f0[1]],

                        pVals[f0[0]],
                        pCoords[f0[0]],
                        snappedPoint[f0[0]],

                        pVals[f0[2]],
                        pCoords[f0[2]],
                        snappedPoint[f0[2]],

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

                        pVals[f0[0]],
                        pCoords[f0[0]],
                        snappedPoint[f0[0]],

                        pVals[f0[1]],
                        pCoords[f0[1]],
                        snappedPoint[f0[1]],

                        pVals[f0[2]],
                        pCoords[f0[2]],
                        snappedPoint[f0[2]],

                        pVals[oppositeI],
                        pCoords[oppositeI],
                        snappedPoint[oppositeI],

                        triPoints
                    );
                }
            }
            else
            {
                forAll(cFaces, cFacei)
                {
                    label facei = cFaces[cFacei];
                    const face& f = mesh_.faces()[facei];

                    label fp0 = mesh_.tetBasePtIs()[facei];

                    // Skip undefined tets
                    if (fp0 < 0)
                    {
                        fp0 = 0;
                        countNotFoundTets++;
                    }

                    label fp = f.fcIndex(fp0);
                    for (label i = 2; i < f.size(); i++)
                    {
                        label nextFp = f.fcIndex(fp);
                        triFace tri(f[fp0], f[fp], f[nextFp]);

                        // Start off from positive volume tet to make sure we
                        // generate outwards pointing tets
                        if (mesh_.faceOwner()[facei] == celli)
                        {
                            generateTriPoints
                            (
                                snappedPoints,

                                pVals[tri[1]],
                                pCoords[tri[1]],
                                snappedPoint[tri[1]],

                                pVals[tri[0]],
                                pCoords[tri[0]],
                                snappedPoint[tri[0]],

                                pVals[tri[2]],
                                pCoords[tri[2]],
                                snappedPoint[tri[2]],

                                cVals[celli],
                                cCoords[celli],
                                snappedCc[celli],

                                triPoints
                            );
                        }
                        else
                        {
                            generateTriPoints
                            (
                                snappedPoints,

                                pVals[tri[0]],
                                pCoords[tri[0]],
                                snappedPoint[tri[0]],

                                pVals[tri[1]],
                                pCoords[tri[1]],
                                snappedPoint[tri[1]],

                                pVals[tri[2]],
                                pCoords[tri[2]],
                                snappedPoint[tri[2]],

                                cVals[celli],
                                cCoords[celli],
                                snappedCc[celli],

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
                triMeshCells.append(celli);
            }
        }
    }

    if (countNotFoundTets > 0)
    {
        WarningIn("Foam::isoSurfaceCell::generateTriPoints(..)")
            << "Could not find " << countNotFoundTets
            << " tet base points, which may lead to inverted triangles."
            << endl;
    }

    triPoints.shrink();
    triMeshCells.shrink();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::isoSurfaceCell::interpolate
(
    const Field<Type>& cCoords,
    const Field<Type>& pCoords
) const
{
    DynamicList<Type> triPoints(3*nCutCells_);
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

    return isoSurface::interpolate
    (
        points().size(),
        triPointMergeMap_,
        triPoints
    );
}


// ************************************************************************* //
