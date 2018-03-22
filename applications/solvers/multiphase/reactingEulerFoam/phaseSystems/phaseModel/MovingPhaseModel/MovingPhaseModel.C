/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "MovingPhaseModel.H"
#include "phaseSystem.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::phi(const volVectorField& U) const
{
    word phiName(IOobject::groupName("phi", this->name()));

    IOobject phiHeader
    (
        phiName,
        U.mesh().time().timeName(),
        U.mesh(),
        IOobject::NO_READ
    );

    if (phiHeader.typeHeaderOk<surfaceScalarField>(true))
    {
        Info<< "Reading face flux field " << phiName << endl;

        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    U.mesh().time().timeName(),
                    U.mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                U.mesh()
            )
        );
    }
    else
    {
        Info<< "Calculating face flux field " << phiName << endl;

        wordList phiTypes
        (
            U.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U.boundaryField(), i)
        {
            if
            (
                isA<fixedValueFvPatchVectorField>(U.boundaryField()[i])
             || isA<slipFvPatchVectorField>(U.boundaryField()[i])
             || isA<partialSlipFvPatchVectorField>(U.boundaryField()[i])
            )
            {
                phiTypes[i] = fixedValueFvPatchScalarField::typeName;
            }
        }

        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    U.mesh().time().timeName(),
                    U.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::flux(U),
                phiTypes
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MovingPhaseModel<BasePhaseModel>::MovingPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, index),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
    ),
    phi_(phi(U_)),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0)
    ),
    alphaRhoPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPhi", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0)
    ),
    DUDt_(nullptr),
    DUDtf_(nullptr),
    divU_(nullptr),
    turbulence_
    (
        phaseCompressibleTurbulenceModel::New
        (
            *this,
            this->thermo().rho(),
            U_,
            alphaRhoPhi_,
            phi_,
            *this
        )
    ),
    continuityErrorFlow_
    (
        IOobject
        (
            IOobject::groupName("continuityErrorFlow", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimDensity/dimTime, 0)
    ),
    continuityErrorSources_
    (
        IOobject
        (
            IOobject::groupName("continuityErrorSources", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimDensity/dimTime, 0)
    ),
    K_(nullptr)
{
    phi_.writeOpt() = IOobject::AUTO_WRITE;

    correctKinematics();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MovingPhaseModel<BasePhaseModel>::~MovingPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correct()
{
    BasePhaseModel::correct();

    this->fluid().MRF().correctBoundaryVelocity(U_);

    volScalarField& rho = this->thermoRef().rho();

    continuityErrorFlow_ = fvc::ddt(*this, rho) + fvc::div(alphaRhoPhi_);

    continuityErrorSources_ = - (this->fluid().fvOptions()(*this, rho)&rho);
}


template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correctKinematics()
{
    BasePhaseModel::correctKinematics();

    if (DUDt_.valid())
    {
        DUDt_.clear();
        DUDt();
    }

    if (DUDtf_.valid())
    {
        DUDtf_.clear();
        DUDtf();
    }

    if (K_.valid())
    {
        K_.clear();
        K();
    }
}


template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correctTurbulence()
{
    BasePhaseModel::correctTurbulence();

    turbulence_->correct();
}


template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correctEnergyTransport()
{
    BasePhaseModel::correctEnergyTransport();

    turbulence_->correctEnergyTransport();
}


template<class BasePhaseModel>
bool Foam::MovingPhaseModel<BasePhaseModel>::stationary() const
{
    return false;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::MovingPhaseModel<BasePhaseModel>::UEqn()
{
    const volScalarField& alpha = *this;
    const volScalarField& rho = this->thermo().rho();

    return
    (
        fvm::ddt(alpha, rho, U_)
      + fvm::div(alphaRhoPhi_, U_)
      + fvm::SuSp(- this->continuityError(), U_)
      + this->fluid().MRF().DDt(alpha*rho, U_)
      + turbulence_->divDevRhoReff(U_)
    );
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::MovingPhaseModel<BasePhaseModel>::UfEqn()
{
    // As the "normal" U-eqn but without the ddt terms

    const volScalarField& alpha = *this;
    const volScalarField& rho = this->thermo().rho();

    return
    (
        fvm::div(alphaRhoPhi_, U_)
      - fvm::Sp(fvc::div(alphaRhoPhi_), U_)
      + fvm::SuSp(- this->continuityErrorSources(), U_)
      + this->fluid().MRF().DDt(alpha*rho, U_)
      + turbulence_->divDevRhoReff(U_)
    );
}


template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::MovingPhaseModel<BasePhaseModel>::U() const
{
    return U_;
}


template<class BasePhaseModel>
Foam::volVectorField&
Foam::MovingPhaseModel<BasePhaseModel>::URef()
{
    return U_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::phi() const
{
    return phi_;
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::MovingPhaseModel<BasePhaseModel>::phiRef()
{
    return phi_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::alphaPhi() const
{
    return alphaPhi_;
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::MovingPhaseModel<BasePhaseModel>::alphaPhiRef()
{
    return alphaPhi_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::alphaRhoPhi() const
{
    return alphaRhoPhi_;
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::MovingPhaseModel<BasePhaseModel>::alphaRhoPhiRef()
{
    return alphaRhoPhi_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::MovingPhaseModel<BasePhaseModel>::DUDt() const
{
    if (!DUDt_.valid())
    {
        DUDt_ = fvc::ddt(U_) + fvc::div(phi_, U_) - fvc::div(phi_)*U_;
    }

    return tmp<volVectorField>(DUDt_());
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::DUDtf() const
{
    if (!DUDtf_.valid())
    {
        DUDtf_ = byDt(phi_ - phi_.oldTime());
    }

    return tmp<surfaceScalarField>(DUDtf_());
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::continuityError() const
{
    return continuityErrorFlow_ + continuityErrorSources_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::continuityErrorFlow() const
{
    return continuityErrorFlow_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::continuityErrorSources() const
{
    return continuityErrorSources_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::K() const
{
    if (!K_.valid())
    {
        K_ =
            new volScalarField
            (
                IOobject::groupName("K", this->name()),
                0.5*magSqr(this->U())
            );
    }

    return tmp<volScalarField>(K_());
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::divU() const
{
    return divU_.valid() ? tmp<volScalarField>(divU_()) : tmp<volScalarField>();
}


template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::divU(tmp<volScalarField> divU)
{
    divU_ = divU;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::mut() const
{
    return turbulence_->mut();
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::muEff() const
{
    return turbulence_->muEff();
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::nut() const
{
    return turbulence_->nut();
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::nuEff() const
{
    return turbulence_->nuEff();
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::kappaEff() const
{
    return turbulence_->kappaEff();
}


template<class BasePhaseModel>
Foam::tmp<Foam::scalarField>
Foam::MovingPhaseModel<BasePhaseModel>::kappaEff(const label patchi) const
{
    return turbulence_->kappaEff(patchi);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::alphaEff() const
{
    return turbulence_->alphaEff();
}


template<class BasePhaseModel>
Foam::tmp<Foam::scalarField>
Foam::MovingPhaseModel<BasePhaseModel>::alphaEff(const label patchi) const
{
    return turbulence_->alphaEff(patchi);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::k() const
{
    return turbulence_->k();
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::pPrime() const
{
    return turbulence_->pPrime();
}


// ************************************************************************* //
