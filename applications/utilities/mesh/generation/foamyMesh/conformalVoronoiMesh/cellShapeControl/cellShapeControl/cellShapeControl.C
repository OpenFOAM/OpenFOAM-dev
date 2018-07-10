/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "cellShapeControl.H"
#include "pointField.H"
#include "scalarField.H"
#include "triadField.H"
#include "cellSizeAndAlignmentControl.H"
#include "searchableSurfaceControl.H"
#include "cellSizeFunction.H"
#include "indexedVertexOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellShapeControl, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellShapeControl::cellShapeControl
(
    const Time& runTime,
    const cvControls& foamyHexMeshControls,
    const searchableSurfaces& allGeometry,
    const conformationSurfaces& geometryToConformTo
)
:
    dictionary
    (
        foamyHexMeshControls.foamyHexMeshDict().subDict("motionControl")
    ),
    runTime_(runTime),
    allGeometry_(allGeometry),
    geometryToConformTo_(geometryToConformTo),
    defaultCellSize_(foamyHexMeshControls.defaultCellSize()),
    minimumCellSize_(foamyHexMeshControls.minimumCellSize()),
    shapeControlMesh_(runTime),
    aspectRatio_(*this),
    sizeAndAlignment_
    (
        runTime,
        subDict("shapeControlFunctions"),
        geometryToConformTo_,
        defaultCellSize_
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellShapeControl::~cellShapeControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalarField Foam::cellShapeControl::cellSize
(
    const pointField& pts
) const
{
    scalarField cellSizes(pts.size());

    forAll(pts, i)
    {
        cellSizes[i] = cellSize(pts[i]);
    }

    return cellSizes;
}


Foam::scalar Foam::cellShapeControl::cellSize(const point& pt) const
{
    barycentric bary;
    cellShapeControlMesh::Cell_handle ch;

    shapeControlMesh_.barycentricCoords(pt, bary, ch);

    scalar size = 0;

    if (shapeControlMesh_.dimension() < 3)
    {
        size = sizeAndAlignment_.cellSize(pt);
    }
    else if (shapeControlMesh_.is_infinite(ch))
    {
//        if (nFarPoints != 0)
//        {
//            for (label pI = 0; pI < 4; ++pI)
//            {
//                if (!ch->vertex(pI)->farPoint())
//                {
//                    size = ch->vertex(pI)->targetCellSize();
//                    return size;
//                }
//            }
//        }

//        cellShapeControlMesh::Vertex_handle nearV =
//            shapeControlMesh_.nearest_vertex_in_cell
//            (
//                toPoint<cellShapeControlMesh::Point>(pt),
//                ch
//            );
//
//        size = nearV->targetCellSize();

        // Find nearest surface. This can be quite slow if there are a lot of
        // surfaces
        size = sizeAndAlignment_.cellSize(pt);
    }
    else
    {
        label nFarPoints = 0;
        for (label pI = 0; pI < 4; ++pI)
        {
            if (ch->vertex(pI)->farPoint())
            {
                nFarPoints++;
            }
        }

        if (nFarPoints != 0)
        {
            for (label pI = 0; pI < 4; ++pI)
            {
                if (!CGAL::indexedVertexOps::uninitialised(ch->vertex(pI)))
                {
                    size = ch->vertex(pI)->targetCellSize();
                    return size;
                }
            }
        }
        else
        {
            forAll(bary, pI)
            {
                size += bary[pI]*ch->vertex(pI)->targetCellSize();
            }
        }
    }

    return size;
}


Foam::tensor Foam::cellShapeControl::cellAlignment(const point& pt) const
{
    barycentric bary;
    cellShapeControlMesh::Cell_handle ch;

    shapeControlMesh_.barycentricCoords(pt, bary, ch);

    tensor alignment = Zero;

    if (shapeControlMesh_.dimension() < 3 || shapeControlMesh_.is_infinite(ch))
    {
        alignment = tensor::I;
    }
    else
    {
        label nFarPoints = 0;
        for (label pI = 0; pI < 4; ++pI)
        {
            if (ch->vertex(pI)->farPoint())
            {
                nFarPoints++;
            }
        }

//        if (nFarPoints != 0)
//        {
//            for (label pI = 0; pI < 4; ++pI)
//            {
//                if (!ch->vertex(pI)->farPoint())
//                {
//                    alignment = ch->vertex(pI)->alignment();
//                }
//            }
//        }
//        else
        {
            triad tri;

            for (label pI = 0; pI < 4; ++pI)
            {
                if (bary[pI] > small)
                {
                    tri += triad(bary[pI]*ch->vertex(pI)->alignment());
                }
            }

            tri.normalize();
            tri.orthogonalize();
            tri = tri.sortxyz();

            alignment = tri;
        }

//        cellShapeControlMesh::Vertex_handle nearV =
//            shapeControlMesh_.nearest_vertex_in_cell
//            (
//                toPoint<cellShapeControlMesh::Point>(pt),
//                ch
//            );
//
//        alignment = nearV->alignment();
    }

    return alignment;
}


void Foam::cellShapeControl::cellSizeAndAlignment
(
    const point& pt,
    scalar& size,
    tensor& alignment
) const
{
    barycentric bary;
    cellShapeControlMesh::Cell_handle ch;

    shapeControlMesh_.barycentricCoords(pt, bary, ch);

    alignment = Zero;
    size = 0;

    if (shapeControlMesh_.dimension() < 3 || shapeControlMesh_.is_infinite(ch))
    {
        // Find nearest surface
        size = sizeAndAlignment_.cellSize(pt);
        alignment = tensor::I;
    }
    else
    {
        label nFarPoints = 0;
        for (label pI = 0; pI < 4; ++pI)
        {
            if (ch->vertex(pI)->farPoint())
            {
                nFarPoints++;
            }
        }

        if (nFarPoints != 0)
        {
            for (label pI = 0; pI < 4; ++pI)
            {
                if (!CGAL::indexedVertexOps::uninitialised(ch->vertex(pI)))
                {
                    size = ch->vertex(pI)->targetCellSize();
                    alignment = ch->vertex(pI)->alignment();
                }
            }
        }
        else
        {
            triad tri;

            for (label pI = 0; pI < 4; ++pI)
            {
                size += bary[pI]*ch->vertex(pI)->targetCellSize();

                if (bary[pI] > small)
                {
                    tri += triad(bary[pI]*ch->vertex(pI)->alignment());
                }
            }

            tri.normalize();
            tri.orthogonalize();
            tri = tri.sortxyz();

            alignment = tri;

//            cellShapeControlMesh::Vertex_handle nearV =
//                shapeControlMesh_.nearest_vertex
//                (
//                    toPoint<cellShapeControlMesh::Point>(pt)
//                );
//
//            alignment = nearV->alignment();
        }
    }

    for (label dir = 0; dir < 3; dir++)
    {
        triad v = alignment;

        if (!v.set(dir) || size == 0)
        {
            // Force orthogonalization of triad.

            scalar dotProd = great;
            if (dir == 0)
            {
                dotProd = v[1] & v[2];

                v[dir] = v[1] ^ v[2];
            }
            if (dir == 1)
            {
                dotProd = v[0] & v[2];

                v[dir] = v[0] ^ v[2];
            }
            if (dir == 2)
            {
                dotProd = v[0] & v[1];

                v[dir] = v[0] ^ v[1];
            }

            v.normalize();
            v.orthogonalize();

            Pout<< "Dot prod = " << dotProd << endl;
            Pout<< "Alignment = " << v << endl;

            alignment = v;

//            FatalErrorInFunction
//                << "Point has bad alignment! "
//                << pt << " " << size << " " << alignment << nl
//                << "Bary Coords = " << bary <<  nl
//                << ch->vertex(0)->info() << nl
//                << ch->vertex(1)->info() << nl
//                << ch->vertex(2)->info() << nl
//                << ch->vertex(3)->info()
//                << abort(FatalError);
        }
    }
}


// ************************************************************************* //
