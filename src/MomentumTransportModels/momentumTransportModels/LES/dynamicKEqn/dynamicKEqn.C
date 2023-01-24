/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "dynamicKEqn.H"
#include "fvModels.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
Foam::tmp<Foam::volScalarField>
dynamicKEqn<BasicMomentumTransportModel>::KK() const
{
    return max
    (
        0.5*(filter_(magSqr(this->U_)) - magSqr(filter_(this->U_))),
        dimensionedScalar(sqr(dimVelocity), small)
    );
}


template<class BasicMomentumTransportModel>
volScalarField dynamicKEqn<BasicMomentumTransportModel>::Ck
(
    const volSymmTensorField& D,
    const volScalarField& KK
) const
{
    const volSymmTensorField LL
    (
        simpleFilter_(dev(filter_(sqr(this->U_)) - (sqr(filter_(this->U_)))))
    );

    const volSymmTensorField MM
    (
        simpleFilter_(-2.0*this->delta()*sqrt(KK)*filter_(D))
    );

    const volScalarField Ck
    (
        simpleFilter_(0.5*(LL && MM))
       /(
            simpleFilter_(magSqr(MM))
          + dimensionedScalar(sqr(MM.dimensions()), vSmall)
        )
    );

    tmp<volScalarField> tfld = 0.5*(mag(Ck) + Ck);
    return tfld();
}


template<class BasicMomentumTransportModel>
volScalarField dynamicKEqn<BasicMomentumTransportModel>::Ce
(
    const volSymmTensorField& D,
    const volScalarField& KK
) const
{
    const volScalarField Ce
    (
        simpleFilter_(this->nuEff()*(filter_(magSqr(D)) - magSqr(filter_(D))))
       /simpleFilter_(pow(KK, 1.5)/(2.0*this->delta()))
    );

    tmp<volScalarField> tfld = 0.5*(mag(Ce) + Ce);
    return tfld();
}


template<class BasicMomentumTransportModel>
volScalarField dynamicKEqn<BasicMomentumTransportModel>::Ce() const
{
    const volSymmTensorField D(dev(symm(fvc::grad(this->U_))));
    return Ce(D, KK());
}


template<class BasicMomentumTransportModel>
void dynamicKEqn<BasicMomentumTransportModel>::correctNut
(
    const volSymmTensorField& D,
    const volScalarField& KK
)
{
    this->nut_ = Ck(D, KK)*sqrt(k_)*this->delta();
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
void dynamicKEqn<BasicMomentumTransportModel>::correctNut()
{
    correctNut(symm(fvc::grad(this->U_)), KK());
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> dynamicKEqn<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
dynamicKEqn<BasicMomentumTransportModel>::dynamicKEqn
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

    simpleFilter_(this->mesh_),
    filterPtr_(LESfilter::New(this->mesh_, this->coeffDict())),
    filter_(filterPtr_())
{
    bound(k_, this->kMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool dynamicKEqn<BasicMomentumTransportModel>::read()
{
    if (LESeddyViscosity<BasicMomentumTransportModel>::read())
    {
        filter_.read(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> dynamicKEqn<BasicMomentumTransportModel>::epsilon() const
{
    return volScalarField::New
    (
        this->groupName("epsilon"),
        Ce()*k()*sqrt(k())/this->delta()
    );
}


template<class BasicMomentumTransportModel>
void dynamicKEqn<BasicMomentumTransportModel>::correct()
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

    LESeddyViscosity<BasicMomentumTransportModel>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volSymmTensorField D(dev(symm(tgradU())));
    const volScalarField G(this->GName(), 2.0*nut*(tgradU() && D));
    tgradU.clear();

    const volScalarField KK(this->KK());

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
    ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(Ce(D, KK)*alpha*rho*sqrt(k_)/this->delta(), k_)
      + kSource()
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);

    correctNut(D, KK);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
