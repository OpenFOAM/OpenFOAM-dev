/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2024 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "solidification.H"
#include "geometricOneField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(solidification, 0);
        addToRunTimeSelectionTable(porosityModel, solidification, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::solidification::solidification
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& coeffDict,
    const word& cellZoneName
)
:
    porosityModel(name, mesh, dict, coeffDict, cellZoneName),
    TName_(coeffDict.lookupOrDefault<word>("T", "T")),
    alphaName_(coeffDict.lookupOrDefault<word>("alpha", "none")),
    rhoName_(coeffDict.lookupOrDefault<word>("rho", "rho")),
    D_(Function1<scalar>::New("D", dimTemperature, dimless/dimTime, coeffDict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::solidification::~solidification()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::solidification::calcTransformModelData()
{}


void Foam::porosityModels::solidification::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    scalarField Udiag(U.size(), 0.0);
    const scalarField& V = mesh_.V();

    apply(Udiag, V, rho, U);

    force = Udiag*U;
}


void Foam::porosityModels::solidification::correct
(
    fvVectorMatrix& UEqn
) const
{
    const volVectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>
        (
            IOobject::groupName(rhoName_, U.group())
        );

        apply(Udiag, V, rho, U);
    }
    else
    {
        apply(Udiag, V, geometricOneField(), U);
    }
}


void Foam::porosityModels::solidification::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const volVectorField& U = UEqn.psi();

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>
        (
            IOobject::groupName(rhoName_, U.group())
        );

        apply(AU, rho, U);
    }
    else
    {
        apply(AU, geometricOneField(), U);
    }
}


// ************************************************************************* //
