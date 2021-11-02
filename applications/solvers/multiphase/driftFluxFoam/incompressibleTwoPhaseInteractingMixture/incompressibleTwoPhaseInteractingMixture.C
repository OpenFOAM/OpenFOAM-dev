/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2021 OpenFOAM Foundation
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

#include "incompressibleTwoPhaseInteractingMixture.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleTwoPhaseInteractingMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleTwoPhaseInteractingMixture::
incompressibleTwoPhaseInteractingMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    twoPhaseMixture(U.mesh()),

    muModel_(mixtureViscosityModel::New(U.mesh(), phase1Name())),
    nucModel_(viscosityModel::New(U.mesh(), phase2Name())),

    rhod_("rho", dimDensity, muModel_()),
    rhoc_("rho", dimDensity, nucModel_()),
    dd_
    (
        "d",
        dimLength,
        muModel_->lookupOrDefault("d", 0.0)
    ),
    alphaMax_(muModel_->lookupOrDefault("alphaMax", 1.0)),

    U_(U),

    mu_
    (
        IOobject
        (
            "mu",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::incompressibleTwoPhaseInteractingMixture::read()
{
    if (twoPhaseMixture::read())
    {
        if (muModel_->read() || nucModel_->read())
        {
            muModel_->lookup("rho") >> rhod_;
            nucModel_->lookup("rho") >> rhoc_;

            dd_ = dimensionedScalar
            (
                "d",
                dimLength,
                muModel_->lookupOrDefault("d", 0)
            );

            alphaMax_ = muModel_->lookupOrDefault( "alphaMax", 1.0);

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
