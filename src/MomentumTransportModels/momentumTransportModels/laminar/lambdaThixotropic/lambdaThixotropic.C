/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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

#include "lambdaThixotropic.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
lambdaThixotropic<BasicMomentumTransportModel>::lambdaThixotropic
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
:
    linearViscousStress<laminarModel<BasicMomentumTransportModel>>
    (
        typeName,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    a_("a", dimless/dimTime, this->coeffDict_),
    b_("b", dimless, this->coeffDict_),
    d_("d", dimless, this->coeffDict_),
    c_("c", pow(dimTime, d_.value() - scalar(1)), this->coeffDict_),
    nu0_("nu0", dimViscosity, this->coeffDict_),
    nuInf_("nuInf", dimViscosity, this->coeffDict_),
    K_(1 - sqrt(nuInf_/nu0_)),

    lambda_
    (
        IOobject
        (
            IOobject::groupName
            (
                typedName("lambda"),
                alphaRhoPhi.group()
            ),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    nu_
    (
        IOobject
        (
            IOobject::groupName
            (
                typedName("nu"),
                alphaRhoPhi.group()
            ),
            this->runTime_.name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField>
lambdaThixotropic<BasicMomentumTransportModel>::calcNu() const
{
    return nuInf_/(sqr(1 - K_*lambda_) + rootVSmall);
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
lambdaThixotropic<BasicMomentumTransportModel>::strainRate() const
{
    return sqrt(2.0)*mag(symm(fvc::grad(this->U())()()));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool lambdaThixotropic<BasicMomentumTransportModel>::read()
{
    if (laminarModel<BasicMomentumTransportModel>::read())
    {
        a_.read(this->coeffDict());
        b_.read(this->coeffDict());
        d_.read(this->coeffDict());

        c_ = dimensionedScalar
        (
            "c",
            pow(dimTime, d_.value() - scalar(1)),
            this->coeffDict_
        );

        nu0_.read(this->coeffDict());
        nuInf_.read(this->coeffDict());

        K_ = (1 - sqrt(nuInf_/nu0_));

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
lambdaThixotropic<BasicMomentumTransportModel>::nuEff() const
{
    return volScalarField::New
    (
        IOobject::groupName("nuEff", this->alphaRhoPhi_.group()),
        nu_
    );
}


template<class BasicMomentumTransportModel>
tmp<scalarField>
lambdaThixotropic<BasicMomentumTransportModel>::nuEff
(
    const label patchi
) const
{
    return nu_.boundaryField()[patchi];
}


template<class BasicMomentumTransportModel>
void lambdaThixotropic<BasicMomentumTransportModel>::correct()
{
    // Local references
    const surfaceScalarField& phi = this->phi_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    tmp<fvScalarMatrix> lambdaEqn
    (
        fvm::ddt(lambda_) + fvm::div(phi, lambda_)
      - fvm::Sp(fvc::div(phi), lambda_)
     ==
        a_*pow(1 - lambda_(), b_)
      - fvm::Sp(c_*pow(strainRate(), d_), lambda_)
      + fvModels.source(lambda_)
    );

    lambdaEqn.ref().relax();
    fvConstraints.constrain(lambdaEqn.ref());
    solve(lambdaEqn);
    fvConstraints.constrain(lambda_);

    lambda_.maxMin(scalar(0), scalar(1));

    nu_ = calcNu();

    laminarModel<BasicMomentumTransportModel>::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarModels
} // End namespace Foam

// ************************************************************************* //
