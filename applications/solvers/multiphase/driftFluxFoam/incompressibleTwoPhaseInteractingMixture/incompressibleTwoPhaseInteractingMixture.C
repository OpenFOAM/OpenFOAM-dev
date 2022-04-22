/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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
#include "relativeVelocityModel.H"
#include "fvcDiv.H"
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
    const surfaceScalarField& phi,
    const uniformDimensionedVectorField& g
)
:
    twoPhaseMixture(U.mesh()),

    U_(U),

    muModel_(mixtureViscosityModel::New(*this)),
    nucModel_(viscosityModel::New(U.mesh(), phase2Name())),

    rhod_("rho", dimDensity, muModel_()),
    rhoc_("rho", dimDensity, nucModel_()),
    dd_
    (
        "d",
        dimLength,
        muModel_->lookupOrDefault("d", 0.0)
    ),
    alphaMax_(lookupOrDefault("alphaMax", 1.0)),

    g_(g),

    MRF_(U.mesh()),

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
    ),

    UdmModel_(relativeVelocityModel::New(*this, *this, g))
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::incompressibleTwoPhaseInteractingMixture::
~incompressibleTwoPhaseInteractingMixture()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volScalarField&
Foam::incompressibleTwoPhaseInteractingMixture::alphad() const
{
    return alpha1();
}


const Foam::volScalarField&
Foam::incompressibleTwoPhaseInteractingMixture::alphac() const
{
    return alpha2();
}


const Foam::mixtureViscosityModel&
Foam::incompressibleTwoPhaseInteractingMixture::muModel() const
{
    return muModel_();
}


const Foam::viscosityModel&
Foam::incompressibleTwoPhaseInteractingMixture::nucModel() const
{
    return nucModel_();
}


const Foam::dimensionedScalar&
Foam::incompressibleTwoPhaseInteractingMixture::rhod() const
{
    return rhod_;
}


const Foam::dimensionedScalar&
Foam::incompressibleTwoPhaseInteractingMixture::rhoc() const
{
    return rhoc_;
};


const Foam::dimensionedScalar&
Foam::incompressibleTwoPhaseInteractingMixture::dd() const
{
    return dd_;
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


const Foam::IOMRFZoneList&
Foam::incompressibleTwoPhaseInteractingMixture::MRF() const
{
    return MRF_;
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseInteractingMixture::mu() const
{
    return mu_;
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleTwoPhaseInteractingMixture::mu(const label patchi) const
{
    return mu_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseInteractingMixture::rho() const
{
    return alpha1()*rhod_ + alpha2()*rhoc_;
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleTwoPhaseInteractingMixture::rho(const label patchi) const
{
    return
        alpha1().boundaryField()[patchi]*rhod_.value()
      + alpha2().boundaryField()[patchi]*rhoc_.value();
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseInteractingMixture::nu() const
{
    return mu_/rho();
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleTwoPhaseInteractingMixture::nu(const label patchi) const
{
    return mu_.boundaryField()[patchi]/rho(patchi);
}


const Foam::volVectorField&
Foam::incompressibleTwoPhaseInteractingMixture::Udm() const
{
    return UdmModel_->Udm();
}


Foam::tmp<Foam::volVectorField>
Foam::incompressibleTwoPhaseInteractingMixture::divTauDm() const
{
    return fvc::div(UdmModel_->tauDm());
}


void Foam::incompressibleTwoPhaseInteractingMixture::correct()
{
    MRF_.correctBoundaryVelocity(U_);
    mu_ = muModel_->mu(rhoc_*nucModel_->nu(), U_);
    UdmModel_->correct();
}


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
