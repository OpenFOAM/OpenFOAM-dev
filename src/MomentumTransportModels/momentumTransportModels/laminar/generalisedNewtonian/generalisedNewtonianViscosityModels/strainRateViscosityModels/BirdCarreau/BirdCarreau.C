/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2024 OpenFOAM Foundation
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

#include "BirdCarreau.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalisedNewtonianViscosityModels
{
    defineTypeNameAndDebug(BirdCarreau, 0);
    addToRunTimeSelectionTable
    (
        generalisedNewtonianViscosityModel,
        BirdCarreau,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarModels::generalisedNewtonianViscosityModels::BirdCarreau::
BirdCarreau
(
    const dictionary& viscosityProperties,
    const Foam::viscosity& viscosity,
    const volVectorField& U
)
:
    strainRateViscosityModel(viscosityProperties, viscosity, U),
    nuInf_("nuInf", dimKinematicViscosity, 0),
    k_("k", dimTime, 0),
    tauStar_( "tauStar", dimKinematicViscosity/dimTime, 0),
    n_("n", dimless, 0),
    a_("a", dimless, 2)
{
    read(viscosityProperties);
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::laminarModels::generalisedNewtonianViscosityModels::BirdCarreau::
read
(
    const dictionary& viscosityProperties
)
{
    strainRateViscosityModel::read(viscosityProperties);

    const dictionary& coeffs =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    nuInf_.read(coeffs);;

    if (coeffs.found("tauStar"))
    {
        tauStar_.read(coeffs);
    }
    else
    {
        k_.read(coeffs);
    }

    n_.read(coeffs);
    a_.readIfPresent(coeffs);

    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarModels::generalisedNewtonianViscosityModels::BirdCarreau::
nu
(
    const volScalarField& nu0,
    const volScalarField& strainRate
) const
{
    return
        nuInf_
      + (nu0 - nuInf_)
       *pow
        (
            scalar(1)
          + pow
            (
                tauStar_.value() > 0
              ? nu0*strainRate/tauStar_
              : k_*strainRate,
                a_
            ),
            (n_ - 1.0)/a_
        );
}


// ************************************************************************* //
