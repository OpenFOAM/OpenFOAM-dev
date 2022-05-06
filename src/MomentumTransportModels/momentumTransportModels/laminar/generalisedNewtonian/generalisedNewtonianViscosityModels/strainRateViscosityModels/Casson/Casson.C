/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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

#include "Casson.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalisedNewtonianViscosityModels
{
    defineTypeNameAndDebug(Casson, 0);
    addToRunTimeSelectionTable
    (
        generalisedNewtonianViscosityModel,
        Casson,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarModels::generalisedNewtonianViscosityModels::Casson::Casson
(
    const dictionary& viscosityProperties,
    const Foam::viscosity& viscosity,
    const volVectorField& U
)
:
    strainRateViscosityModel(viscosityProperties, viscosity, U),
    m_("m", dimViscosity, 0),
    tau0_("tau0", dimViscosity/dimTime, 0),
    nuMin_("nuMin", dimViscosity, 0),
    nuMax_("nuMax", dimViscosity, 0)
{
    read(viscosityProperties);
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::laminarModels::generalisedNewtonianViscosityModels::Casson::read
(
    const dictionary& viscosityProperties
)
{
    strainRateViscosityModel::read(viscosityProperties);

    const dictionary& coeffs =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    m_.read(coeffs);
    tau0_.read(coeffs);
    nuMin_.read(coeffs);
    nuMax_.read(coeffs);

    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarModels::generalisedNewtonianViscosityModels::Casson::
nu
(
    const volScalarField& nu0,
    const volScalarField& strainRate
) const
{
    return max
    (
        nuMin_,
        min
        (
            nuMax_,
            sqr
            (
                sqrt
                (
                    tau0_
                   /max
                    (
                        strainRate,
                        dimensionedScalar(dimless/dimTime, vSmall)
                    )
                ) + sqrt(m_)
            )
        )
    );
}


// ************************************************************************* //
