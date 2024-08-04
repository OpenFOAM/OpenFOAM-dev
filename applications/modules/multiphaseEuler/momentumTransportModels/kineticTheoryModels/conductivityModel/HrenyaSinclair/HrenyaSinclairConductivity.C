/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "HrenyaSinclairConductivity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace conductivityModels
{
    defineTypeNameAndDebug(HrenyaSinclair, 0);

    addToRunTimeSelectionTable
    (
        conductivityModel,
        HrenyaSinclair,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::kineticTheoryModels::conductivityModels::HrenyaSinclair::readCoeffs
(
    const dictionary& coeffDict
)
{
    L_.readIfPresent(coeffDict);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::conductivityModels::HrenyaSinclair::HrenyaSinclair
(
    const dictionary& coeffDict
)
:
    conductivityModel(coeffDict),
    L_("L", dimensionSet(0, 1, 0, 0, 0), coeffDict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::conductivityModels::HrenyaSinclair::
~HrenyaSinclair()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::conductivityModels::HrenyaSinclair::kappa
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

    volScalarField lambda
    (
        scalar(1) + da/(6*sqrt(2.0)*(alpha1 + 1e-5))/L_
    );

    return rho1*da*sqrt(Theta)*
    (
        2*sqr(alpha1)*g0*(1 + e)/sqrtPi
      + (9.0/8.0)*sqrtPi*g0*0.25*sqr(1 + e)*(2*e - 1)*sqr(alpha1)
       /(49.0/16.0 - 33*e/16.0)
      + (15.0/16.0)*sqrtPi*alpha1*(0.5*sqr(e) + 0.25*e - 0.75 + lambda)
       /((49.0/16.0 - 33*e/16.0)*lambda)
      + (25.0/64.0)*sqrtPi
       /((1 + e)*(49.0/16.0 - 33*e/16.0)*lambda*g0)
    );
}


// ************************************************************************* //
