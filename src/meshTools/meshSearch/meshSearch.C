/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "meshSearch.H"
#include "meshSearchBoundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshSearch, 0);
}


// * * * * * * * * * * * * * * Protected Constructors  * * * * * * * * * * * //

Foam::meshSearch::meshSearch(const polyMesh& mesh)
:
    DemandDrivenMeshObject
    <
        polyMesh,
        DeletableMeshObject,
        meshSearch
    >(mesh),
    cellTree_
    (
       treeDataCell(false, mesh),
       meshSearchBoundBox::New(mesh).bb(),
       8,
       10,
       scalar(5)
    )
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

const Foam::meshSearch& Foam::meshSearch::New
(
    const polyMesh& mesh,
    const pointInCellShapes cellShapes
)
{
    if
    (
        cellShapes == pointInCellShapes::faceDiagonalTris
     || cellShapes == pointInCellShapes::tets
    )
    {
        // Force construction of face diagonals
        (void)mesh.tetBasePtIs();
    }

    return
        DemandDrivenMeshObject
        <
            polyMesh,
            DeletableMeshObject,
            meshSearch
        >::New(mesh);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshSearch::~meshSearch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::meshSearch::findNearestCell(const point& p) const
{
    pointIndexHit info =
        cellTree().findNearest
        (
            p,
            magSqr(cellTree().bb().max() - cellTree().bb().min())
        );

    if (!info.hit())
    {
        info = cellTree().findNearest(p, Foam::sqr(great));
    }

    return info.index();
}


Foam::label Foam::meshSearch::findNearestCellNoTree
(
    const polyMesh& mesh,
    const point& p
)
{
    const pointField& cellCentres = mesh.cellCentres();

    label nearestCelli = -1;
    scalar nearestDistSqr = vGreat;

    forAll(cellCentres, celli)
    {
        const scalar distSqr = magSqr(cellCentres[celli] - p);

        if (distSqr < nearestDistSqr)
        {
            nearestCelli = celli;
            nearestDistSqr = distSqr;
        }
    }

    return nearestCelli;
}


Foam::label Foam::meshSearch::findNearestFace
(
    const point& p
) const
{
    // Find the nearest cell
    const label celli = findNearestCell(p);
    if (celli < 0) return -1;

    const pointField& faceCentres = mesh().faceCentres();

    label nearestFacei = -1;
    scalar nearestDistSqr = vGreat;

    // Check all faces of the nearest cell
    const cell& cFaces = mesh().cells()[celli];

    forAll(cFaces, cFacei)
    {
        const label facei = cFaces[cFacei];

        const scalar distSqr = magSqr(faceCentres[facei] - p);

        if (distSqr < nearestDistSqr)
        {
            nearestFacei = facei;
            nearestDistSqr = distSqr;
        }
    }

    return nearestFacei;
}


Foam::label Foam::meshSearch::findCell
(
    const point& p,
    const pointInCellShapes cellShapes
) const
{
    return cellTree().findInside(p, cellShapes);
}


Foam::label Foam::meshSearch::findCellNoTree
(
    const polyMesh& mesh,
    const point& p,
    const pointInCellShapes cellShapes
)
{
    // Find the nearest cell
    const label celli = findNearestCellNoTree(mesh, p);
    if (celli < 0) return -1;

    // If point is in the nearest cell return
    if (pointInCell(p, mesh, celli, cellShapes))
    {
        return celli;
    }

    // If the point is not in the nearest cell then search all cells
    for (label celli = 0; celli < mesh.nCells(); celli++)
    {
        if (pointInCell(p, mesh, celli, cellShapes))
        {
            return celli;
        }
    }

    return -1;
}


// ************************************************************************* //
