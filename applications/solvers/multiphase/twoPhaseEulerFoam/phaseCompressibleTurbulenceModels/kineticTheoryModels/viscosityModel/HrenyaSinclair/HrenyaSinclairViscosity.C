/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "HrenyaSinclairViscosity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace viscosityModels
{
    defineTypeNameAndDebug(HrenyaSinclair, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        HrenyaSinclair,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::viscosityModels::HrenyaSinclair::HrenyaSinclair
(
    const dictionary& dict
)
:
    viscosityModel(dict),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
    L_("L", dimensionSet(0, 1, 0, 0, 0), coeffDict_.lookup("L"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::viscosityModels::HrenyaSinclair::~HrenyaSinclair()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::viscosityModels::HrenyaSinclair::nu
(
    const volScalarField& alpha1,
    const volScalarField& Theta,
    const volScalarField& g0,
    const volScalarField& rho1,
    const volScalarField& da,
    const dimensionedScalar& e
) const
{
    const scalar sqrtPi = sqrt(constant::mathematical::pi);

    volScalarField lamda
    (
        scalar(1) + da/(6.0*sqrt(2.0)*(alpha1 + scalar(1.0e-5)))/L_
    );

    return da*sqrt(Theta)*
    (
        (4.0/5.0)*sqr(alpha1)*g0*(1.0 + e)/sqrtPi
      + (1.0/15.0)*sqrtPi*g0*(1.0 + e)*(3.0*e - 1)*sqr(alpha1)/(3.0-e)
      + (1.0/6.0)*sqrtPi*alpha1*(0.5*lamda + 0.25*(3.0*e - 1.0))
       /(0.5*(3.0 - e)*lamda)
      + (10/96.0)*sqrtPi/((1.0 + e)*0.5*(3.0 - e)*g0*lamda)
    );
}


bool Foam::kineticTheoryModels::viscosityModels::HrenyaSinclair::read()
{
    coeffDict_ <<= dict_.subDict(typeName + "Coeffs");

    L_.readIfPresent(coeffDict_);

    return true;
}


// ************************************************************************* //
