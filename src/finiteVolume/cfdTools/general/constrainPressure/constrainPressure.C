/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "constrainPressure.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "geometricOneField.H"
#include "fixedFluxPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class RhoType, class RAUType, class MRFType>
void Foam::constrainPressure
(
    volScalarField& p,
    const RhoType& rho,
    const volVectorField& U,
    const surfaceScalarField& phiHbyA,
    const RAUType& rhorAU,
    const MRFType& MRF
)
{
    const fvMesh& mesh = p.mesh();

    volScalarField::Boundary& pBf = p.boundaryFieldRef();

    const volVectorField::Boundary& UBf = U.boundaryField();
    const surfaceScalarField::Boundary& phiHbyABf =
        phiHbyA.boundaryField();
    const typename RAUType::Boundary& rhorAUBf =
        rhorAU.boundaryField();
    const surfaceVectorField::Boundary& SfBf =
        mesh.Sf().boundaryField();
    const surfaceScalarField::Boundary& magSfBf =
        mesh.magSf().boundaryField();

    forAll(pBf, patchi)
    {
        if (isA<fixedFluxPressureFvPatchScalarField>(pBf[patchi]))
        {
            refCast<fixedFluxPressureFvPatchScalarField>
            (
                pBf[patchi]
            ).updateCoeffs
            (
                (
                    phiHbyABf[patchi]
                  - rho.boundaryField()[patchi]
                   *MRF.relative(SfBf[patchi] & UBf[patchi], patchi)
                )
               /(magSfBf[patchi]*rhorAUBf[patchi])
            );
        }
    }
}


template<class RAUType>
void Foam::constrainPressure
(
    volScalarField& p,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phiHbyA,
    const RAUType& rAU
)
{
    constrainPressure(p, rho, U, phiHbyA, rAU, NullMRF());
}


template<class RAUType, class MRFType>
void Foam::constrainPressure
(
    volScalarField& p,
    const volVectorField& U,
    const surfaceScalarField& phiHbyA,
    const RAUType& rAU,
    const MRFType& MRF
)
{
    constrainPressure(p, geometricOneField(), U, phiHbyA, rAU, MRF);
}


template<class RAUType>
void Foam::constrainPressure
(
    volScalarField& p,
    const volVectorField& U,
    const surfaceScalarField& phiHbyA,
    const RAUType& rAU
)
{
    constrainPressure(p, U, phiHbyA, rAU, NullMRF());
}


// ************************************************************************* //
