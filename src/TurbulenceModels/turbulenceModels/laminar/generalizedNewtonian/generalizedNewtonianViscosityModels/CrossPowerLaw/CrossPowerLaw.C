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

#include "CrossPowerLaw.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalizedNewtonianViscosityModels
{
    defineTypeNameAndDebug(CrossPowerLaw, 0);

    addToRunTimeSelectionTable
    (
        generalizedNewtonianViscosityModel,
        CrossPowerLaw,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarModels::generalizedNewtonianViscosityModels::CrossPowerLaw::
CrossPowerLaw
(
    const dictionary& viscosityProperties
)
:
    generalizedNewtonianViscosityModel(viscosityProperties),
    nuInf_("nuInf", dimViscosity, 0),
    m_("m", dimTime, 0),
    tauStar_( "tauStar", dimViscosity/dimTime, 0),
    n_("n", dimless, 0)
{
    read(viscosityProperties);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::laminarModels::generalizedNewtonianViscosityModels::CrossPowerLaw::
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
        coeffs.lookup("m") >> m_;
    }

    coeffs.lookup("n") >> n_;

    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarModels::generalizedNewtonianViscosityModels::CrossPowerLaw::
nu
(
    const volScalarField& nu0,
    const volScalarField& strainRate
) const
{
    return
       nuInf_
     + (nu0 - nuInf_)
       /(
           scalar(1)
         + pow
           (
                tauStar_.value() > 0
              ? nu0*strainRate/tauStar_
              : m_*strainRate,
                n_
            )
        );
}


// ************************************************************************* //
