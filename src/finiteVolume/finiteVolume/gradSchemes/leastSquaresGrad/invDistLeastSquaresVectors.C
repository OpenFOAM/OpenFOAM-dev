/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "leastSquaresVectors.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(leastSquaresVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::leastSquaresVectors::leastSquaresVectors(const fvMesh& mesh)
:
    MeshObject<fvMesh, Foam::MoveableMeshObject, leastSquaresVectors>(mesh),
    pVectors_
    (
        IOobject
        (
            "LeastSquaresP",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    ),
    nVectors_
    (
        IOobject
        (
            "LeastSquaresN",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    )
{
    calcLeastSquaresVectors();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::leastSquaresVectors::~leastSquaresVectors()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::leastSquaresVectors::calcLeastSquaresVectors()
{
    if (debug)
    {
        Info<< "leastSquaresVectors::calcLeastSquaresVectors() :"
            << "Calculating least square gradient vectors"
            << endl;
    }

    const fvMesh& mesh = mesh_;

    // Set local references to mesh data
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const volVectorField& C = mesh.C();

    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh_.nCells(), symmTensor::zero);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        vector d = C[nei] - C[own];
        symmTensor wdd = sqr(d)/magSqr(d);
        dd[own] += wdd;
        dd[nei] += wdd;
    }


    surfaceVectorField::GeometricBoundaryField& blsP =
        pVectors_.boundaryField();

    forAll(blsP, patchi)
    {
        const fvsPatchVectorField& patchLsP = blsP[patchi];

        const fvPatch& p = patchLsP.patch();
        const labelUList& faceCells = p.patch().faceCells();

        // Build the d-vectors
        vectorField pd(p.delta());

        forAll(pd, patchFacei)
        {
            const vector& d = pd[patchFacei];

            dd[faceCells[patchFacei]] += sqr(d)/magSqr(d);
        }
    }


    // Invert the dd tensor
    const symmTensorField invDd(inv(dd));


    // Revisit all faces and calculate the pVectors_ and nVectors_ vectors
    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        vector d = C[nei] - C[own];

        pVectors_[facei] = (invDd[own] & d)/magSqr(d);
        nVectors_[facei] = -(invDd[nei] & d)/magSqr(d);
    }

    forAll(blsP, patchi)
    {
        fvsPatchVectorField& patchLsP = blsP[patchi];

        const fvPatch& p = patchLsP.patch();
        const labelUList& faceCells = p.faceCells();

        // Build the d-vectors
        vectorField pd(p.delta());

        forAll(pd, patchFacei)
        {
            const vector& d = pd[patchFacei];

            patchLsP[patchFacei] = (invDd[faceCells[patchFacei]] & d)/magSqr(d);
        }
    }

    if (debug)
    {
        Info<< "leastSquaresVectors::calcLeastSquaresVectors() :"
            << "Finished calculating least square gradient vectors"
            << endl;
    }
}


bool Foam::leastSquaresVectors::movePoints()
{
    calcLeastSquaresVectors();
    return true;
}


// ************************************************************************* //
