/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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

#include "v2f.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField> v2f<BasicMomentumTransportModel>::boundEpsilon()
{
    tmp<volScalarField> tCmuk2(CmuKEps_*sqr(k_));
    epsilon_ = max(epsilon_, tCmuk2()/(this->nutMaxCoeff_*this->nu()));
    return tCmuk2;
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> v2f<BasicMomentumTransportModel>::Ts() const
{
    return max(k_/epsilon_, 6.0*sqrt(this->nu()/epsilon_));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> v2f<BasicMomentumTransportModel>::Ls() const
{
    return
        CL_*max(pow(k_, 1.5)
       /epsilon_, Ceta_*pow025(pow3(this->nu())/epsilon_));
}


template<class BasicMomentumTransportModel>
void v2f<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = min(boundEpsilon()/epsilon_, this->Cmu_*v2_*Ts());
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
v2f<BasicMomentumTransportModel>::v2f
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
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),
    v2fBase(),

    Cmu_("Cmu", this->coeffDict(), 0.22),
    CmuKEps_("CmuKEps", this->coeffDict(), 0.09),
    C1_("C1", this->coeffDict(), 1.4),
    C2_("C2", this->coeffDict(), 0.3),
    CL_("CL", this->coeffDict(), 0.23),
    Ceta_("Ceta", this->coeffDict(), 70.0),
    Ceps2_("Ceps2", this->coeffDict(), 1.9),
    Ceps3_("Ceps3", this->coeffDict(), -0.33),
    sigmaK_("sigmaK", this->coeffDict(), 1.0),
    sigmaEps_("sigmaEps", this->coeffDict(), 1.3),

    k_
    (
        IOobject
        (
            this->groupName("k"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            this->groupName("epsilon"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    v2_
    (
        IOobject
        (
            this->groupName("v2"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    f_
    (
        IOobject
        (
            this->groupName("f"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    v2Min_(dimensionedScalar(v2_.dimensions(), small)),
    fMin_(dimensionedScalar(f_.dimensions(), 0))
{
    bound(k_, this->kMin_);
    boundEpsilon();
    bound(v2_, v2Min_);
    bound(f_, fMin_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool v2f<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        CmuKEps_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        CL_.readIfPresent(this->coeffDict());
        Ceta_.readIfPresent(this->coeffDict());
        Ceps2_.readIfPresent(this->coeffDict());
        Ceps3_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void v2f<BasicMomentumTransportModel>::correct()
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
    volScalarField& nut = this->nut_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    // Use N=6 so that f=0 at walls
    const dimensionedScalar N("N", dimless, 6.0);

    const volTensorField gradU(fvc::grad(U));
    const volScalarField S2(2*magSqr(dev(symm(gradU))));

    const volScalarField G(this->GName(), nut*S2);
    const volScalarField Ts(this->Ts());
    const volScalarField L2(typedName("L2"), sqr(Ls()));
    const volScalarField v2fAlpha
    (
        typedName("alpha"),
        1.0/Ts*((C1_ - N)*v2_ - 2.0/3.0*k_*(C1_ - 1.0))
    );

    const volScalarField Ceps1
    (
        typedName("Ceps1"),
        1.4*(1.0 + 0.05*min(sqrt(k_/v2_), scalar(100.0)))
    );

    // Update epsilon (and possibly G) at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        Ceps1*alpha*rho*G/Ts
      - fvm::SuSp(((2.0/3.0)*Ceps1 + Ceps3_)*alpha*rho*divU, epsilon_)
      - fvm::Sp(Ceps2_*alpha*rho/Ts, epsilon_)
      + fvModels.source(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvConstraints.constrain(epsilon_);
    boundEpsilon();


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(alpha*rho*epsilon_/k_, k_)
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);


    // Relaxation function equation
    tmp<fvScalarMatrix> fEqn
    (
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/L2, f_)
      - 1.0/L2/k_*(v2fAlpha - C2_*G)
    );

    fEqn.ref().relax();
    fvConstraints.constrain(fEqn.ref());
    solve(fEqn);
    fvConstraints.constrain(f_);
    bound(f_, fMin_);


    // Turbulence stress normal to streamlines equation
    tmp<fvScalarMatrix> v2Eqn
    (
        fvm::ddt(alpha, rho, v2_)
      + fvm::div(alphaRhoPhi, v2_)
      - fvm::laplacian(alpha*rho*DkEff(), v2_)
      ==
        alpha*rho*min(k_*f_, C2_*G - v2fAlpha)
      - fvm::Sp(N*alpha*rho*epsilon_/k_, v2_)
      + fvModels.source(alpha, rho, v2_)
    );

    v2Eqn.ref().relax();
    fvConstraints.constrain(v2Eqn.ref());
    solve(v2Eqn);
    fvConstraints.constrain(v2_);
    bound(v2_, v2Min_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
