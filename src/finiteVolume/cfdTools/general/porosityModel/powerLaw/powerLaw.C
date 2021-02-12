/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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
#include "powerLaw.H"
#include "geometricOneField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(powerLaw, 0);
        addToRunTimeSelectionTable(porosityModel, powerLaw, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::powerLaw::powerLaw
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porosityModel(name, modelType, mesh, dict, cellZoneName),
    C0_(coeffs_.lookup<scalar>("C0")),
    C1_(coeffs_.lookup<scalar>("C1")),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "rho"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::powerLaw::~powerLaw()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::powerLaw::calcTransformModelData()
{
    // nothing to be transformed
}


void Foam::porosityModels::powerLaw::calcForce
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


void Foam::porosityModels::powerLaw::correct
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


void Foam::porosityModels::powerLaw::correct
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


bool Foam::porosityModels::powerLaw::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
