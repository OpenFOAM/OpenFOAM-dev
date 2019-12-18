/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "consumptionSpeed.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(consumptionSpeed, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::consumptionSpeed::consumptionSpeed
(
    const dictionary& dict
)
:   omega0_(dict.lookup<scalar>("omega0")),
    eta_(dict.lookup<scalar>("eta")),
    sigmaExt_(dict.lookup<scalar>("sigmaExt")),
    omegaMin_(dict.lookup<scalar>("omegaMin"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::consumptionSpeed::~consumptionSpeed()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::consumptionSpeed::omega0Sigma
(
    scalar sigma,
    scalar a
) const
{
    scalar omega0 = 0.0;

    if (sigma < sigmaExt_)
    {
        omega0 = max
        (
            a*omega0_*(1.0 - exp(eta_*max(sigma, 0.0))),
            omegaMin_
        ) ;
    }

    return omega0;
}


Foam::tmp<Foam::volScalarField> Foam::consumptionSpeed::omega0Sigma
(
    const volScalarField& sigma
)
{
    tmp<volScalarField> tomega0
    (
        volScalarField::New
        (
            "omega0",
            sigma.mesh(),
            dimensionedScalar
            (
                dimensionSet(1, -2, -1, 0, 0, 0, 0),
                0
            )
        )
    );

    volScalarField& omega0 = tomega0.ref();

    volScalarField::Internal& iomega0 = omega0;

    forAll(iomega0, celli)
    {
        iomega0[celli] = omega0Sigma(sigma[celli], 1.0);
    }

    volScalarField::Boundary& bomega0 = omega0.boundaryFieldRef();

    forAll(bomega0, patchi)
    {
        forAll(bomega0[patchi], facei)
        {
            bomega0[patchi][facei] =
                omega0Sigma
                (
                    sigma.boundaryField()[patchi][facei],
                    1.0
                );
        }
    }

    return tomega0;
}


void  Foam::consumptionSpeed::read(const dictionary& dict)
{
    dict.lookup("omega0") >> omega0_ ;
    dict.lookup("eta") >> eta_ ;
    dict.lookup("sigmaExt") >> sigmaExt_;
    dict.lookup("omegaMin") >> omegaMin_;
}


// ************************************************************************* //
