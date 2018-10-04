/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalizedNewtonianViscosityModels
{
    defineTypeNameAndDebug(BirdCarreau, 0);
    addToRunTimeSelectionTable
    (
        generalizedNewtonianViscosityModel,
        BirdCarreau,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarModels::generalizedNewtonianViscosityModels::BirdCarreau::
BirdCarreau
(
    const dictionary& viscosityProperties
)
:
    generalizedNewtonianViscosityModel(viscosityProperties),
    nuInf_("nuInf", dimViscosity, 0),
    k_("k", dimTime, 0),
    tauStar_( "tauStar", dimViscosity/dimTime, 0),
    n_("n", dimless, 0),
    a_("a", dimless, 2)
{
    read(viscosityProperties);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::laminarModels::generalizedNewtonianViscosityModels::BirdCarreau::
read
(
    const dictionary& viscosityProperties
)
{
    generalizedNewtonianViscosityModel::read(viscosityProperties);

    const dictionary& coeffs =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    coeffs.lookup("nuInf") >> nuInf_;

    if (coeffs.found("tauStar"))
    {
        coeffs.lookup("tauStar") >> tauStar_;
    }
    else
    {
        coeffs.lookup("k") >> k_;
    }

    coeffs.lookup("n") >> n_;

    a_ = coeffs.lookupOrDefault
    (
        "a",
        dimensionedScalar("a", dimless, 2)
    );

    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarModels::generalizedNewtonianViscosityModels::BirdCarreau::
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
