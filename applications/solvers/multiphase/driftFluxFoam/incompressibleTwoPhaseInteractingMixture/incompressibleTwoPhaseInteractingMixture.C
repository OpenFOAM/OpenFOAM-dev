/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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
#include "mixtureViscosityModel.H"
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
    volVectorField& U,
    const surfaceScalarField& phi
)
:
    twoPhaseMixture(U.mesh()),

    U_(U),

    nucModel_(viscosityModel::New(U.mesh(), phase2Name())),
    muModel_(mixtureViscosityModel::New(*this)),

    rhoc_("rho", dimDensity, nucModel_()),
    rhod_("rho", dimDensity, muModel_()),

    alphaMax_(lookupOrDefault("alphaMax", 1.0)),

    rho_
    (
        IOobject
        (
            "rho",
            U_.time().name(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("rho", dimDensity, 0)
    ),

    nu_
    (
        IOobject
        (
            "nu",
            U_.time().name(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar(dimViscosity, 0),
        calculatedFvPatchScalarField::typeName
    )
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::incompressibleTwoPhaseInteractingMixture::
~incompressibleTwoPhaseInteractingMixture()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volScalarField&
Foam::incompressibleTwoPhaseInteractingMixture::alphac() const
{
    return alpha2();
}


const Foam::volScalarField&
Foam::incompressibleTwoPhaseInteractingMixture::alphad() const
{
    return alpha1();
}


const Foam::dimensionedScalar&
Foam::incompressibleTwoPhaseInteractingMixture::rhoc() const
{
    return rhoc_;
};


const Foam::dimensionedScalar&
Foam::incompressibleTwoPhaseInteractingMixture::rhod() const
{
    return rhod_;
}


Foam::scalar
Foam::incompressibleTwoPhaseInteractingMixture::alphaMax() const
{
    return alphaMax_;
}


const Foam::volVectorField&
Foam::incompressibleTwoPhaseInteractingMixture::U() const
{
    return U_;
}

const Foam::volScalarField&
Foam::incompressibleTwoPhaseInteractingMixture::rho() const
{
    return rho_;
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseInteractingMixture::nu() const
{
    return nu_;
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleTwoPhaseInteractingMixture::nu(const label patchi) const
{
    return nu_.boundaryField()[patchi];
}


void Foam::incompressibleTwoPhaseInteractingMixture::correct()
{
    rho_ = alpha1()*rhod_ + alpha2()*rhoc_;
    nu_ = muModel_->mu(rhoc_*nucModel_->nu(), U_)/rho_;
}


bool Foam::incompressibleTwoPhaseInteractingMixture::read()
{
    if (twoPhaseMixture::read())
    {
        if (muModel_->read() || nucModel_->read())
        {
            nucModel_->lookup("rho") >> rhoc_;
            muModel_->lookup("rho") >> rhod_;

            alphaMax_ = muModel_->lookupOrDefault("alphaMax", 1.0);

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
