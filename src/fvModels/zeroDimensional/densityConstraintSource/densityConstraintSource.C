/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "densityConstraintSource.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
namespace zeroDimensional
{
    defineTypeNameAndDebug(densityConstraintSource, 0);
    addToRunTimeSelectionTable(fvModel, densityConstraintSource, dictionary);
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::zeroDimensional::densityConstraintSource::readCoeffs()
{
    rho_.reset(Function1<scalar>::New("rho", coeffs()).ptr());
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::zeroDimensional::densityConstraintSource::dmdt() const
{
    const fluidThermo& thermo =
        mesh().lookupObject<fluidThermo>(physicalProperties::typeName);

    const scalar t = mesh().time().userTimeValue();
    const dimensionedScalar& deltaT = mesh().time().deltaT();

    const dimensionedScalar rhoTarget(dimPressure, rho_->value(t));

    const volScalarField::Internal& rho = thermo.rho();

    return (rhoTarget - rho)/deltaT;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::zeroDimensional::densityConstraintSource::densityConstraintSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    constraintSource(name, modelType, mesh, dict),
    rho_(nullptr)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::zeroDimensional::densityConstraintSource::~densityConstraintSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::zeroDimensional::densityConstraintSource::read
(
    const dictionary& dict
)
{
    if (fvModel::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
