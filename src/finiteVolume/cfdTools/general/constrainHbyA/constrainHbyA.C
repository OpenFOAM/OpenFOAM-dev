/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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

#include "constrainHbyA.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fixedFluxExtrapolatedPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::constrainHbyA
(
    const tmp<volVectorField>& tHbyA,
    const volVectorField& U,
    const volScalarField& p
)
{
    tmp<volVectorField> tHbyANew;

    if (tHbyA.isTmp())
    {
        tHbyANew = tHbyA;
        tHbyANew.ref().rename(IOobject::groupName("HbyA", U.group()));
    }
    else
    {
        tHbyANew = volVectorField::New
        (
            IOobject::groupName("HbyA", U.group()),
            tHbyA
        );
    }

    volVectorField& HbyA = tHbyANew.ref();
    volVectorField::Boundary& HbyAbf = HbyA.boundaryFieldRef();

    forAll(U.boundaryField(), patchi)
    {
        if
        (
           !U.boundaryField()[patchi].assignable()
        && !isA<fixedFluxExtrapolatedPressureFvPatchScalarField>
            (
                p.boundaryField()[patchi]
            )
        )
        {
            HbyAbf[patchi] = U.boundaryField()[patchi];
        }
    }

    return tHbyANew;
}


Foam::tmp<Foam::surfaceScalarField> Foam::constrainPhiHbyA
(
    const tmp<surfaceScalarField>& tphiHbyA,
    const volVectorField& U,
    const volScalarField& p
)
{
    tmp<surfaceScalarField> tphiHbyANew;

    if (tphiHbyA.isTmp())
    {
        tphiHbyANew = tphiHbyA;
        tphiHbyANew.ref().rename(IOobject::groupName("phiHbyA", U.group()));
    }
    else
    {
        tphiHbyANew = surfaceScalarField::New
        (
            IOobject::groupName("phiHbyA", U.group()),
            tphiHbyA
        );
    }

    surfaceScalarField& phiHbyA = tphiHbyANew.ref();
    surfaceScalarField::Boundary& phiHbyAbf = phiHbyA.boundaryFieldRef();

    forAll(U.boundaryField(), patchi)
    {
        if
        (
           !U.boundaryField()[patchi].assignable()
        && !isA<fixedFluxExtrapolatedPressureFvPatchScalarField>
            (
                p.boundaryField()[patchi]
            )
        )
        {
            phiHbyAbf[patchi] =
                U.mesh().Sf().boundaryField()[patchi]
              & U.boundaryField()[patchi];
        }
    }

    return tphiHbyANew;
}


Foam::tmp<Foam::surfaceScalarField> Foam::constrainPhid
(
    const tmp<surfaceScalarField>& tphid,
    const volScalarField& p
)
{
    surfaceScalarField& phid = tphid.ref();
    surfaceScalarField::Boundary& phidBf = phid.boundaryFieldRef();

    const volScalarField::Boundary& pBf = p.boundaryField();

    forAll(pBf, patchi)
    {
        if (isA<fixedFluxPressureFvPatchScalarField>(pBf[patchi]))
        {
            phidBf[patchi] = 0;
        }
    }

    return tphid;
}


// ************************************************************************* //
