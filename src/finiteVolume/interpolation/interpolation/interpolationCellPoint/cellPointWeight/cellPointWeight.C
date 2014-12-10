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

#include "cellPointWeight.H"
#include "polyMesh.H"
#include "polyMeshTetDecomposition.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::cellPointWeight::debug(debug::debugSwitch("cellPointWeight", 0));

Foam::scalar Foam::cellPointWeight::tol(SMALL);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::cellPointWeight::findTetrahedron
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI
)
{
    if (debug)
    {
        Pout<< nl << "Foam::cellPointWeight::findTetrahedron" << nl
            << "position = " << position << nl
            << "cellI = " << cellI << endl;
    }

    List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
    (
        mesh,
        cellI
    );

    const faceList& pFaces = mesh.faces();
    const scalar cellVolume = mesh.cellVolumes()[cellI];

    forAll(cellTets, tetI)
    {
        const tetIndices& tetIs = cellTets[tetI];

        const face& f = pFaces[tetIs.face()];

        // Barycentric coordinates of the position
        scalar det = tetIs.tet(mesh).barycentric(position, weights_);

        if (mag(det/cellVolume) > tol)
        {
            const scalar& u = weights_[0];
            const scalar& v = weights_[1];
            const scalar& w = weights_[2];

            if
            (
                (u + tol > 0)
             && (v + tol > 0)
             && (w + tol > 0)
             && (u + v + w < 1 + tol)
            )
            {
                faceVertices_[0] = f[tetIs.faceBasePt()];
                faceVertices_[1] = f[tetIs.facePtA()];
                faceVertices_[2] = f[tetIs.facePtB()];

                return;
            }
        }
    }

    // A suitable point in a tetrahedron was not found, find the
    // nearest.

    scalar minNearDist = VGREAT;

    label nearestTetI = -1;

    forAll(cellTets, tetI)
    {
        const tetIndices& tetIs = cellTets[tetI];

        scalar nearDist = tetIs.tet(mesh).nearestPoint(position).distance();

        if (nearDist < minNearDist)
        {
            minNearDist = nearDist;

            nearestTetI = tetI;
        }
    }

    if (debug)
    {
        Pout<< "cellPointWeight::findTetrahedron" << nl
            << "    Tetrahedron search failed; using closest tet to point "
            << position << nl
            << "    cell: "
            << cellI << nl
            << endl;
    }


    const tetIndices& tetIs = cellTets[nearestTetI];

    const face& f = pFaces[tetIs.face()];

    // Barycentric coordinates of the position, ignoring if the
    // determinant is suitable.  If not, the return from barycentric
    // to weights_ is safe.
    tetIs.tet(mesh).barycentric(position, weights_);

    faceVertices_[0] = f[tetIs.faceBasePt()];
    faceVertices_[1] = f[tetIs.facePtA()];
    faceVertices_[2] = f[tetIs.facePtB()];
}


void Foam::cellPointWeight::findTriangle
(
    const polyMesh& mesh,
    const vector& position,
    const label faceI
)
{
    if (debug)
    {
        Pout<< "\nbool Foam::cellPointWeight::findTriangle" << nl
            << "position = " << position << nl
            << "faceI = " << faceI << endl;
    }

    List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
    (
        mesh,
        faceI,
        mesh.faceOwner()[faceI]
    );

    const scalar faceAreaSqr = magSqr(mesh.faceAreas()[faceI]);

    const face& f =  mesh.faces()[faceI];

    forAll(faceTets, tetI)
    {
        const tetIndices& tetIs = faceTets[tetI];

        List<scalar> triWeights(3);

        // Barycentric coordinates of the position
        scalar det = tetIs.faceTri(mesh).barycentric(position, triWeights);

        if (0.25*mag(det)/faceAreaSqr > tol)
        {
            const scalar& u = triWeights[0];
            const scalar& v = triWeights[1];

            if
            (
                (u + tol > 0)
             && (v + tol > 0)
             && (u + v < 1 + tol)
            )
            {
                // Weight[0] is for the cell centre.
                weights_[0] = 0;
                weights_[1] = triWeights[0];
                weights_[2] = triWeights[1];
                weights_[3] = triWeights[2];

                faceVertices_[0] = f[tetIs.faceBasePt()];
                faceVertices_[1] = f[tetIs.facePtA()];
                faceVertices_[2] = f[tetIs.facePtB()];

                return;
            }
        }
    }

    // A suitable point in a triangle was not found, find the nearest.

    scalar minNearDist = VGREAT;

    label nearestTetI = -1;

    forAll(faceTets, tetI)
    {
        const tetIndices& tetIs = faceTets[tetI];

        scalar nearDist = tetIs.faceTri(mesh).nearestPoint(position).distance();

        if (nearDist < minNearDist)
        {
            minNearDist = nearDist;

            nearestTetI = tetI;
        }
    }

    if (debug)
    {
        Pout<< "cellPointWeight::findTriangle" << nl
            << "    Triangle search failed; using closest tri to point "
            << position << nl
            << "    face: "
            << faceI << nl
            << endl;
    }

    const tetIndices& tetIs = faceTets[nearestTetI];

    // Barycentric coordinates of the position, ignoring if the
    // determinant is suitable.  If not, the return from barycentric
    // to triWeights is safe.

    List<scalar> triWeights(3);

    tetIs.faceTri(mesh).barycentric(position, triWeights);

    // Weight[0] is for the cell centre.
    weights_[0] = 0;
    weights_[1] = triWeights[0];
    weights_[2] = triWeights[1];
    weights_[3] = triWeights[2];

    faceVertices_[0] = f[tetIs.faceBasePt()];
    faceVertices_[1] = f[tetIs.facePtA()];
    faceVertices_[2] = f[tetIs.facePtB()];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellPointWeight::cellPointWeight
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label faceI
)
:
    cellI_(cellI),
    weights_(4),
    faceVertices_(3)
{
    if (faceI < 0)
    {
        // Face data not supplied
        findTetrahedron(mesh, position, cellI);
    }
    else
    {
        // Face data supplied
        findTriangle(mesh, position, faceI);
    }
}


// ************************************************************************* //
