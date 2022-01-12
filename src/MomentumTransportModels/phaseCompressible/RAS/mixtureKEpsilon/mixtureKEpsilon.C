/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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

#include "mixtureKEpsilon.H"
#include "fvModels.H"
#include "bound.H"
#include "phaseSystem.H"
#include "dispersedDragModel.H"
#include "dispersedVirtualMassModel.H"
#include "fixedValueFvPatchFields.H"
#include "inletOutletFvPatchFields.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
mixtureKEpsilon<BasicMomentumTransportModel>::mixtureKEpsilon
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

    liquidTurbulencePtr_(nullptr),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            C2_.value()
        )
    ),
    Cp_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cp",
            this->coeffDict_,
            0.25
        )
    ),
    alphap_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphap",
            this->coeffDict_,
            1
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
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
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


template<class BasicMomentumTransportModel>
wordList mixtureKEpsilon<BasicMomentumTransportModel>::epsilonBoundaryTypes
(
    const volScalarField& epsilon
) const
{
    const volScalarField::Boundary& ebf = epsilon.boundaryField();

    wordList ebt = ebf.types();

    forAll(ebf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(ebf[patchi]))
        {
            ebt[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    return ebt;
}


template<class BasicMomentumTransportModel>
void mixtureKEpsilon<BasicMomentumTransportModel>::correctInletOutlet
(
    volScalarField& vsf,
    const volScalarField& refVsf
) const
{
    volScalarField::Boundary& bf = vsf.boundaryFieldRef();
    const volScalarField::Boundary& refBf =
        refVsf.boundaryField();

    forAll(bf, patchi)
    {
        if
        (
            isA<inletOutletFvPatchScalarField>(bf[patchi])
         && isA<inletOutletFvPatchScalarField>(refBf[patchi])
        )
        {
            refCast<inletOutletFvPatchScalarField>
            (bf[patchi]).refValue() =
            refCast<const inletOutletFvPatchScalarField>
            (refBf[patchi]).refValue();
        }
    }
}


template<class BasicMomentumTransportModel>
void mixtureKEpsilon<BasicMomentumTransportModel>::initMixtureFields()
{
    if (rhom_.valid()) return;

    // Local references to gas-phase properties
    const volScalarField& kg = this->k_;
    const volScalarField& epsilong = this->epsilon_;

    // Local references to liquid-phase properties
    mixtureKEpsilon<BasicMomentumTransportModel>& turbc =
        this->liquidTurbulence();
    const volScalarField& kl = turbc.k_;
    const volScalarField& epsilonl = turbc.epsilon_;

    word startTimeName
    (
        this->runTime_.timeName(this->runTime_.startTime().value())
    );

    Ct2_.set
    (
        new volScalarField
        (
            IOobject
            (
                "Ct2",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            Ct2()
        )
    );

    rhom_.set
    (
        new volScalarField
        (
            IOobject
            (
                "rhom",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            rhom()
        )
    );

    km_.set
    (
        new volScalarField
        (
            IOobject
            (
                "km",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mix(kl, kg),
            kl.boundaryField().types()
        )
    );
    correctInletOutlet(km_(), kl);

    epsilonm_.set
    (
        new volScalarField
        (
            IOobject
            (
                "epsilonm",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mix(epsilonl, epsilong),
            epsilonBoundaryTypes(epsilonl)
        )
    );
    correctInletOutlet(epsilonm_(), epsilonl);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool mixtureKEpsilon<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        Cp_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void mixtureKEpsilon<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
mixtureKEpsilon<BasicMomentumTransportModel>&
mixtureKEpsilon<BasicMomentumTransportModel>::liquidTurbulence() const
{
    if (!liquidTurbulencePtr_)
    {
        const volVectorField& U = this->U_;

        const phaseModel& gas = refCast<const phaseModel>(this->properties());
        const phaseSystem& fluid = gas.fluid();
        const phaseModel& liquid = fluid.otherPhase(gas);

        liquidTurbulencePtr_ =
           &const_cast<mixtureKEpsilon<BasicMomentumTransportModel>&>
            (
                U.db().lookupObject
                <
                    mixtureKEpsilon<BasicMomentumTransportModel>
                >
                (
                    IOobject::groupName
                    (
                        momentumTransportModel::typeName,
                        liquid.name()
                    )
                )
            );
    }

    return *liquidTurbulencePtr_;
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> mixtureKEpsilon<BasicMomentumTransportModel>::Ct2() const
{
    const mixtureKEpsilon<BasicMomentumTransportModel>& liquidTurbulence =
        this->liquidTurbulence();

    const phaseModel& gas = refCast<const phaseModel>(this->properties());
    const phaseSystem& fluid = gas.fluid();
    const phaseModel& liquid = fluid.otherPhase(gas);

    const dragModels::dispersedDragModel& drag =
        fluid.lookupInterfacialModel<dragModels::dispersedDragModel>
        (dispersedPhaseInterface(gas, liquid));

    const volScalarField& alphag = this->alpha_;

    volScalarField magUr(mag(liquidTurbulence.U() - this->U()));

    volScalarField beta
    (
        (6*this->Cmu_/(4*sqrt(3.0/2.0)))
       *drag.K()/liquid.rho()
       *(liquidTurbulence.k_/liquidTurbulence.epsilon_)
    );
    volScalarField Ct0((3 + beta)/(1 + beta + 2*gas.rho()/liquid.rho()));
    volScalarField fAlphad((180 + (-4.71e3 + 4.26e4*alphag)*alphag)*alphag);

    return sqr(1 + (Ct0 - 1)*exp(-fAlphad));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
mixtureKEpsilon<BasicMomentumTransportModel>::rholEff() const
{
    const phaseModel& gas = refCast<const phaseModel>(this->properties());
    const phaseSystem& fluid = gas.fluid();
    return fluid.otherPhase(gas).rho();
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
mixtureKEpsilon<BasicMomentumTransportModel>::rhogEff() const
{
    const phaseModel& gas = refCast<const phaseModel>(this->properties());
    const phaseSystem& fluid = gas.fluid();
    const phaseModel& liquid = fluid.otherPhase(gas);

    const virtualMassModels::dispersedVirtualMassModel& virtualMass =
        fluid.lookupInterfacialModel
        <virtualMassModels::dispersedVirtualMassModel>
        (dispersedPhaseInterface(gas, liquid));

    return gas.rho() + virtualMass.Cvm()*fluid.otherPhase(gas).rho();
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> mixtureKEpsilon<BasicMomentumTransportModel>::rhom() const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return alphal*rholEff() + alphag*rhogEff();
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> mixtureKEpsilon<BasicMomentumTransportModel>::mix
(
    const volScalarField& fc,
    const volScalarField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return (alphal*rholEff()*fc + alphag*rhogEff()*fd)/rhom_();
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> mixtureKEpsilon<BasicMomentumTransportModel>::mixU
(
    const volScalarField& fc,
    const volScalarField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return
        (alphal*rholEff()*fc + alphag*rhogEff()*Ct2_()*fd)
       /(alphal*rholEff() + alphag*rhogEff()*Ct2_());
}


template<class BasicMomentumTransportModel>
tmp<surfaceScalarField> mixtureKEpsilon<BasicMomentumTransportModel>::mixFlux
(
    const surfaceScalarField& fc,
    const surfaceScalarField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    surfaceScalarField alphalf(fvc::interpolate(alphal));
    surfaceScalarField alphagf(fvc::interpolate(alphag));

    surfaceScalarField rholEfff(fvc::interpolate(rholEff()));
    surfaceScalarField rhogEfff(fvc::interpolate(rhogEff()));

    return
       (alphalf*rholEfff*fc + alphagf*rhogEfff*fvc::interpolate(Ct2_())*fd)
      /(alphalf*rholEfff + alphagf*rhogEfff*fvc::interpolate(Ct2_()));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
mixtureKEpsilon<BasicMomentumTransportModel>::bubbleG() const
{
    const mixtureKEpsilon<BasicMomentumTransportModel>& liquidTurbulence =
        this->liquidTurbulence();

    const phaseModel& gas = refCast<const phaseModel>(this->properties());
    const phaseSystem& fluid = gas.fluid();
    const phaseModel& liquid = fluid.otherPhase(gas);

    const dragModels::dispersedDragModel& drag =
        fluid.lookupInterfacialModel<dragModels::dispersedDragModel>
        (dispersedPhaseInterface(gas, liquid));

    volScalarField magUr(mag(liquidTurbulence.U() - this->U()));

    // Lahey model
    tmp<volScalarField> bubbleG
    (
        Cp_
       *pos(alphap_ - gas)*liquid*liquid.rho()
       *(
            pow3(magUr)
          + pow(drag.CdRe()*liquid.thermo().nu()/gas.d(), 4.0/3.0)
           *pow(magUr, 5.0/3.0)
        )
       *gas
       /gas.d()
    );

    // Simple model
    // tmp<volScalarField> bubbleG
    // (
    //     Cp_*liquid*drag.K()*sqr(magUr)
    // );

    return bubbleG;
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
mixtureKEpsilon<BasicMomentumTransportModel>::kSource() const
{
    return fvm::Su(bubbleG()/rhom_(), km_());
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
mixtureKEpsilon<BasicMomentumTransportModel>::epsilonSource() const
{
    return fvm::Su(C3_*epsilonm_()*bubbleG()/(rhom_()*km_()), epsilonm_());
}


template<class BasicMomentumTransportModel>
void mixtureKEpsilon<BasicMomentumTransportModel>::correct()
{
    const phaseModel& gas = refCast<const phaseModel>(this->properties());
    const phaseSystem& fluid = gas.fluid();

    // Only solve the mixture turbulence for the gas-phase
    if (&gas != &fluid.phases()[0])
    {
        // This is the liquid phase but check the model for the gas-phase
        // is consistent
        this->liquidTurbulence();

        return;
    }

    if (!this->turbulence_)
    {
        return;
    }

    // Initialise the mixture fields if they have not yet been constructed
    initMixtureFields();

    // Local references to gas-phase properties
    tmp<surfaceScalarField> phig = this->phi();
    const volVectorField& Ug = this->U_;
    const volScalarField& alphag = this->alpha_;
    volScalarField& kg = this->k_;
    volScalarField& epsilong = this->epsilon_;
    volScalarField& nutg = this->nut_;

    // Local references to liquid-phase properties
    mixtureKEpsilon<BasicMomentumTransportModel>& liquidTurbulence =
        this->liquidTurbulence();
    tmp<surfaceScalarField> phil = liquidTurbulence.phi();
    const volVectorField& Ul = liquidTurbulence.U_;
    const volScalarField& alphal = liquidTurbulence.alpha_;
    volScalarField& kl = liquidTurbulence.k_;
    volScalarField& epsilonl = liquidTurbulence.epsilon_;
    volScalarField& nutl = liquidTurbulence.nut_;

    // Local references to mixture properties
    volScalarField& rhom = rhom_();
    volScalarField& km = km_();
    volScalarField& epsilonm = epsilonm_();

    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    // Update the effective mixture density
    rhom = this->rhom();

    // Mixture flux
    const surfaceScalarField phim("phim", mixFlux(phil, phig));

    // Mixture velocity divergence
    const volScalarField divUm
    (
        mixU
        (
            fvc::div(fvc::absolute(phil, Ul)),
            fvc::div(fvc::absolute(phig, Ug))
        )
    );

    tmp<volScalarField> Gc;
    {
        tmp<volTensorField> tgradUl = fvc::grad(Ul);
        Gc = tmp<volScalarField>
        (
            new volScalarField
            (
                this->GName(),
                nutl*(tgradUl() && dev(twoSymm(tgradUl())))
            )
        );
        tgradUl.clear();

        // Update k, epsilon and G at the wall
        kl.boundaryFieldRef().updateCoeffs();
        epsilonl.boundaryFieldRef().updateCoeffs();

        Gc.ref().checkOut();
    }

    tmp<volScalarField> Gd;
    {
        tmp<volTensorField> tgradUg = fvc::grad(Ug);
        Gd = tmp<volScalarField>
        (
            new volScalarField
            (
                this->GName(),
                nutg*(tgradUg() && dev(twoSymm(tgradUg())))
            )
        );
        tgradUg.clear();

        // Update k, epsilon and G at the wall
        kg.boundaryFieldRef().updateCoeffs();
        epsilong.boundaryFieldRef().updateCoeffs();

        Gd.ref().checkOut();
    }

    // Mixture turbulence generation
    const volScalarField Gm(mix(Gc, Gd));

    // Mixture turbulence viscosity
    const volScalarField nutm(mixU(nutl, nutg));

    // Update the mixture k and epsilon boundary conditions
    km == mix(kl, kg);
    bound(km, this->kMin_);
    epsilonm == mix(epsilonl, epsilong);
    bound(epsilonm, this->epsilonMin_);

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilonm)
      + fvm::div(phim, epsilonm)
      - fvm::Sp(fvc::div(phim), epsilonm)
      - fvm::laplacian(DepsilonEff(nutm), epsilonm)
     ==
        C1_*Gm*epsilonm/km
      - fvm::SuSp(((2.0/3.0)*C1_)*divUm, epsilonm)
      - fvm::Sp(C2_*epsilonm/km, epsilonm)
      + epsilonSource()
      + fvModels.source(epsilonm)
    );

    epsEqn.ref().relax();
    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilonm.boundaryFieldRef());
    solve(epsEqn);
    fvConstraints.constrain(epsilonm);
    bound(epsilonm, this->epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kmEqn
    (
        fvm::ddt(km)
      + fvm::div(phim, km)
      - fvm::Sp(fvc::div(phim), km)
      - fvm::laplacian(DkEff(nutm), km)
     ==
        Gm
      - fvm::SuSp((2.0/3.0)*divUm, km)
      - fvm::Sp(epsilonm/km, km)
      + kSource()
      + fvModels.source(km)
    );

    kmEqn.ref().relax();
    fvConstraints.constrain(kmEqn.ref());
    solve(kmEqn);
    fvConstraints.constrain(km);
    bound(km, this->kMin_);
    km.correctBoundaryConditions();

    const volScalarField Cc2(rhom/(alphal*rholEff() + alphag*rhogEff()*Ct2_()));
    kl = Cc2*km;
    kl.correctBoundaryConditions();
    epsilonl = Cc2*epsilonm;
    epsilonl.correctBoundaryConditions();
    liquidTurbulence.correctNut();

    Ct2_() = Ct2();
    kg = Ct2_()*kl;
    kg.correctBoundaryConditions();
    epsilong = Ct2_()*epsilonl;
    epsilong.correctBoundaryConditions();
    nutg = Ct2_()*(liquidTurbulence.nu()/this->nu())*nutl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
