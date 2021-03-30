/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcFlux.H"

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

        forAll(U.boundaryField(), patchi)
        {
            if (!U.boundaryField()[patchi].assignable())
            {
                phiTypes[patchi] = fixedValueFvPatchScalarField::typeName;
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
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index),
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
        dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), 0)
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
        dimensionedScalar(dimensionSet(1, 0, -1, 0, 0), 0)
    ),
    Uf_(nullptr),
    DUDt_(nullptr),
    DUDtf_(nullptr),
    divU_(nullptr),
    turbulence_
    (
        phaseCompressible::momentumTransportModel::New
        (
            *this,
            this->thermo().rho(),
            U_,
            alphaRhoPhi_,
            phi_,
            *this
        )
    ),
    thermophysicalTransport_
    (
        PhaseThermophysicalTransportModel
        <
            phaseCompressible::momentumTransportModel,
            transportThermoModel
        >::New(turbulence_, this->thermo_)
    ),
    continuityError_
    (
        IOobject
        (
            IOobject::groupName("continuityError", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar(dimDensity/dimTime, 0)
    ),
    K_(nullptr)
{
    phi_.writeOpt() = IOobject::AUTO_WRITE;

    if (fluid.mesh().dynamic())
    {
        Uf_ = new surfaceVectorField
        (
            IOobject
            (
                IOobject::groupName("Uf", this->name()),
                fluid.mesh().time().timeName(),
                fluid.mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U_)
        );
    }

    correctKinematics();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MovingPhaseModel<BasePhaseModel>::~MovingPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correctContinuityError
(
    const volScalarField& source
)
{
    volScalarField& rho = this->thermoRef().rho();

    continuityError_ = fvc::ddt(*this, rho) + fvc::div(alphaRhoPhi_) - source;
}


template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correct()
{
    BasePhaseModel::correct();
    this->fluid().MRF().correctBoundaryVelocity(U_);
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
        K_.ref() = 0.5*magSqr(this->U());
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
    thermophysicalTransport_->correct();
}


template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correctUf()
{
    const fvMesh& mesh = this->fluid().mesh();

    if (mesh.dynamic())
    {
        Uf_.ref() = fvc::interpolate(U_);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        Uf_.ref() += n*(fvc::absolute(phi_, U_)/mesh.magSf() - (n & Uf_()));

        surfaceVectorField::Boundary& UfBf = Uf_.ref().boundaryFieldRef();
        const volVectorField::Boundary& UBf = U_.boundaryField();

        forAll(mesh.boundary(), patchi)
        {
            // Remove the flux correction on AMI patches to compensate for
            // AMI non-conservation error
            if (isA<cyclicAMIFvPatch>(mesh.boundary()[patchi]))
            {
                UfBf[patchi] = UBf[patchi];
            }
        }
    }
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
      + fvm::SuSp(-this->continuityError(), U_)
      + this->fluid().MRF().DDt(alpha*rho, U_)
      + turbulence_->divDevTau(U_)
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
      + fvm::SuSp(fvc::ddt(*this, rho) - this->continuityError(), U_)
      + this->fluid().MRF().DDt(alpha*rho, U_)
      + turbulence_->divDevTau(U_)
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
Foam::tmp<Foam::surfaceVectorField>
Foam::MovingPhaseModel<BasePhaseModel>::Uf() const
{
    return
        Uf_.valid()
      ? tmp<surfaceVectorField>(Uf_())
      : tmp<surfaceVectorField>();
}


template<class BasePhaseModel>
Foam::surfaceVectorField&
Foam::MovingPhaseModel<BasePhaseModel>::UfRef()
{
    if (Uf_.valid())
    {
        return Uf_.ref();
    }
    else
    {
        FatalErrorInFunction
            << "Uf has not been allocated."
            << exit(FatalError);

        return const_cast<surfaceVectorField&>(surfaceVectorField::null());
    }
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
        const tmp<surfaceScalarField> taphi(fvc::absolute(phi_, U_));
        const surfaceScalarField& aphi(taphi());
        DUDt_ = fvc::ddt(U_) + fvc::div(aphi, U_) - fvc::div(aphi)*U_;
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
    return continuityError_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::K() const
{
    if (!K_.valid())
    {
        K_ = volScalarField::New
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


template<class BasePhaseModel>
Foam::tmp<Foam::scalarField>
Foam::MovingPhaseModel<BasePhaseModel>::kappaEff(const label patchi) const
{
    return thermophysicalTransport_->kappaEff(patchi);
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MovingPhaseModel<BasePhaseModel>::divq(volScalarField& he) const
{
    return thermophysicalTransport_->divq(he);
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MovingPhaseModel<BasePhaseModel>::divj(volScalarField& Yi) const
{
    return thermophysicalTransport_->divj(Yi);
}


// ************************************************************************* //
