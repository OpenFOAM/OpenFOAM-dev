/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
#include "defineDebugSwitch.H"
#include "polyMesh.H"
#include "polyMeshTetDecomposition.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellPointWeight, 0);
}

Foam::scalar Foam::cellPointWeight::tol(small);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::cellPointWeight::findTetrahedron
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
{
    if (debug)
    {
        Pout<< nl << "Foam::cellPointWeight::findTetrahedron" << nl
            << "position = " << position << nl
            << "celli = " << celli << endl;
    }

    List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
    (
        mesh,
        celli
    );

    const scalar cellVolume = mesh.cellVolumes()[celli];

    forAll(cellTets, tetI)
    {
        const tetIndices& tetIs = cellTets[tetI];

        // Barycentric coordinates of the position
        scalar det = tetIs.tet(mesh).pointToBarycentric(position, weights_);

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

                faceVertices_ = tetIs.faceTriIs(mesh);

                return;
            }
        }
    }

    // A suitable point in a tetrahedron was not found, find the
    // nearest.

    scalar minNearDist = vGreat;

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
            << celli << nl
            << endl;
    }


    const tetIndices& tetIs = cellTets[nearestTetI];

    // Barycentric coordinates of the position, ignoring if the
    // determinant is suitable.  If not, the return from barycentric
    // to weights_ is safe.
    weights_ = tetIs.tet(mesh).pointToBarycentric(position);

    faceVertices_ = tetIs.faceTriIs(mesh);
}


void Foam::cellPointWeight::findTriangle
(
    const polyMesh& mesh,
    const vector& position,
    const label facei
)
{
    if (debug)
    {
        Pout<< "\nbool Foam::cellPointWeight::findTriangle" << nl
            << "position = " << position << nl
            << "facei = " << facei << endl;
    }

    List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
    (
        mesh,
        facei,
        mesh.faceOwner()[facei]
    );

    const scalar faceAreaSqr = sqr(mesh.magFaceAreas()[facei]);

    forAll(faceTets, tetI)
    {
        const tetIndices& tetIs = faceTets[tetI];

        // Barycentric coordinates of the position
        barycentric2D triWeights;
        const scalar det =
            tetIs.faceTri(mesh).pointToBarycentric(position, triWeights);

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

                faceVertices_ = tetIs.faceTriIs(mesh);

                return;
            }
        }
    }

    // A suitable point in a triangle was not found, find the nearest.

    scalar minNearDist = vGreat;

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
            << facei << nl
            << endl;
    }

    const tetIndices& tetIs = faceTets[nearestTetI];

    // Barycentric coordinates of the position, ignoring if the
    // determinant is suitable.  If not, the return from barycentric
    // to triWeights is safe.

    const barycentric2D triWeights =
        tetIs.faceTri(mesh).pointToBarycentric(position);

    // Weight[0] is for the cell centre.
    weights_[0] = 0;
    weights_[1] = triWeights[0];
    weights_[2] = triWeights[1];
    weights_[3] = triWeights[2];

    faceVertices_ = tetIs.faceTriIs(mesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellPointWeight::cellPointWeight
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const label facei
)
:
    celli_(celli)
{
    if (facei < 0)
    {
        // Face data not supplied
        findTetrahedron(mesh, position, celli);
    }
    else
    {
        // Face data supplied
        findTriangle(mesh, position, facei);
    }
}


// ************************************************************************* //
