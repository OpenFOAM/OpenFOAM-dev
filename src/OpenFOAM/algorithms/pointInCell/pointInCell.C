/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "pointInCell.H"
#include "polyMeshTetDecomposition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::pointInCellFacePlanes
(
    const point& p,
    const polyMesh& mesh,
    const label celli
)
{
    const labelList& cFaces = mesh.cells()[celli];

    forAll(cFaces, cFacei)
    {
        const label facei = cFaces[cFacei];
        const bool isOwn = mesh.faceOwner()[facei] == celli;

        const vector normal =
            isOwn ? mesh.faceAreas()[facei] : - mesh.faceAreas()[facei];

        const vector proj = p - mesh.faceCentres()[facei];

        if ((normal & proj) > 0)
        {
            return false;
        }
    }

    return true;
}


bool Foam::pointInCellFaceCentreTris
(
    const point& p,
    const polyMesh& mesh,
    const label celli
)
{
    const cell& cFaces = mesh.cells()[celli];

    forAll(cFaces, cFacei)
    {
        const label facei = cFaces[cFacei];
        const face& f = mesh.faces()[facei];
        const point& fCentre = mesh.faceCentres()[facei];
        const bool isOwn = mesh.faceOwner()[facei] == celli;

        forAll(f, fPointi)
        {
            label pointi;
            label nextPointi;
            if (isOwn)
            {
                pointi = f[fPointi];
                nextPointi = f.nextLabel(fPointi);
            }
            else
            {
                pointi = f.nextLabel(fPointi);
                nextPointi = f[fPointi];
            }

            const triPointRef faceTri
            (
                mesh.points()[pointi],
                mesh.points()[nextPointi],
                fCentre
            );

            const vector proj = p - faceTri.centre();

            if ((faceTri.area() & proj) > 0)
            {
                return false;
            }
        }
    }

    return true;
}


bool Foam::pointInCellFaceDiagTris
(
    const point& p,
    const polyMesh& mesh,
    const label celli
)
{
    const cell& cFaces = mesh.cells()[celli];

    forAll(cFaces, cFacei)
    {
        const label facei = cFaces[cFacei];
        const face& f = mesh.faces()[facei];

        for (label tetPti = 1; tetPti < f.size() - 1; tetPti++)
        {
            const triPointRef faceTri =
                tetIndices(celli, facei, tetPti).faceTri(mesh);

            const vector proj = p - faceTri.centre();

            if ((faceTri.area() & proj) > 0)
            {
                return false;
            }
        }
    }

    return true;
}


bool Foam::pointInCellTets
(
    const point& p,
    const polyMesh& mesh,
    const label celli
)
{
    return polyMeshTetDecomposition::findTet(mesh, celli, p).face() != -1;
}


bool Foam::pointInCell
(
    const point& p,
    const polyMesh& mesh,
    const label celli,
    const pointInCellShapes cellShapes
)
{
    switch (cellShapes)
    {
        case pointInCellShapes::facePlanes:
            return pointInCellFacePlanes(p, mesh, celli);
        case pointInCellShapes::faceCentreTris:
            return pointInCellFaceCentreTris(p, mesh, celli);
        case pointInCellShapes::faceDiagonalTris:
            return pointInCellFaceDiagTris(p, mesh, celli);
        case pointInCellShapes::tets:
            return pointInCellTets(p, mesh, celli);
    }

    return false;

}


// ************************************************************************* //
