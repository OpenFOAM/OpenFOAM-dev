/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "dynamicLagrangian.H"
#include "fvModels.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void dynamicLagrangian<BasicMomentumTransportModel>::correctNut
(
    const tmp<volTensorField>& gradU
)
{
    this->nut_ = (flm_/fmm_)*sqr(this->delta())*mag(dev(symm(gradU)));
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
void dynamicLagrangian<BasicMomentumTransportModel>::correctNut()
{
    correctNut(fvc::grad(this->U_));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
dynamicLagrangian<BasicMomentumTransportModel>::dynamicLagrangian
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    LESeddyViscosity<BasicMomentumTransportModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    flm_
    (
        IOobject
        (
            IOobject::groupName("flm", this->alphaRhoPhi_.group()),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    fmm_
    (
        IOobject
        (
            IOobject::groupName("fmm", this->alphaRhoPhi_.group()),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    theta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "theta",
            this->coeffDict_,
            1.5
        )
    ),

    simpleFilter_(U.mesh()),
    filterPtr_(LESfilter::New(U.mesh(), this->coeffDict())),
    filter_(filterPtr_()),

    flm0_("flm0", flm_.dimensions(), 0.0),
    fmm0_("fmm0", fmm_.dimensions(), vSmall)
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool dynamicLagrangian<BasicMomentumTransportModel>::read()
{
    if (LESeddyViscosity<BasicMomentumTransportModel>::read())
    {
        filter_.read(this->coeffDict());
        theta_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void dynamicLagrangian<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    LESeddyViscosity<BasicMomentumTransportModel>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField S(dev(symm(gradU)));
    volScalarField magS(mag(S));

    volVectorField Uf(filter_(U));
    volSymmTensorField Sf(dev(symm(fvc::grad(Uf))));
    volScalarField magSf(mag(Sf));

    volSymmTensorField L(dev(filter_(sqr(U)) - (sqr(filter_(U)))));
    volSymmTensorField M
    (
        2.0*sqr(this->delta())*(filter_(magS*S) - 4.0*magSf*Sf)
    );

    volScalarField invT
    (
        alpha*rho*(1.0/(theta_.value()*this->delta()))*pow(flm_*fmm_, 1.0/8.0)
    );

    volScalarField LM(L && M);

    fvScalarMatrix flmEqn
    (
        fvm::ddt(alpha, rho, flm_)
      + fvm::div(alphaRhoPhi, flm_)
     ==
        invT*LM
      - fvm::Sp(invT, flm_)
      + fvModels.source(alpha, rho, flm_)
    );

    flmEqn.relax();
    fvConstraints.constrain(flmEqn);
    flmEqn.solve();
    fvConstraints.constrain(flm_);
    bound(flm_, flm0_);

    volScalarField MM(M && M);

    fvScalarMatrix fmmEqn
    (
        fvm::ddt(alpha, rho, fmm_)
      + fvm::div(alphaRhoPhi, fmm_)
     ==
        invT*MM
      - fvm::Sp(invT, fmm_)
      + fvModels.source(alpha, rho, fmm_)
    );

    fmmEqn.relax();
    fvConstraints.constrain(fmmEqn);
    fmmEqn.solve();
    fvConstraints.constrain(fmm_);
    bound(fmm_, fmm0_);

    correctNut(gradU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
