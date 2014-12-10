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

#include "fvcReconstruct.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> reconstructMag(const surfaceScalarField& ssf)
{
    const fvMesh& mesh = ssf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    tmp<volScalarField> treconField
    (
        new volScalarField
        (
            IOobject
            (
                "reconstruct("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "0",
                ssf.dimensions()/dimArea,
                scalar(0)
            ),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    scalarField& rf = treconField();

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        rf[own] += (Sf[facei] & (Cf[facei] - C[own]))*ssf[facei]/magSf[facei];
        rf[nei] -= (Sf[facei] & (Cf[facei] - C[nei]))*ssf[facei]/magSf[facei];
    }

    const surfaceScalarField::GeometricBoundaryField& bsf = ssf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvsPatchScalarField& psf = bsf[patchi];

        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];
        const vectorField& pSf = Sf.boundaryField()[patchi];
        const scalarField& pMagSf = magSf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            label own = pOwner[pFacei];
            rf[own] +=
                (pSf[pFacei] & (pCf[pFacei] - C[own]))
               *psf[pFacei]/pMagSf[pFacei];
        }
    }

    rf /= mesh.V();

    treconField().correctBoundaryConditions();

    return treconField;
}


tmp<volScalarField> reconstructMag(const tmp<surfaceScalarField>& tssf)
{
    tmp<volScalarField> tvf
    (
        fvc::reconstructMag(tssf())
    );
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
