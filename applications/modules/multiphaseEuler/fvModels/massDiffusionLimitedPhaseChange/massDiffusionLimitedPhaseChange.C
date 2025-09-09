/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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

#include "massDiffusionLimitedPhaseChange.H"
#include "fvmSup.H"
#include "multiphaseEuler.H"
#include "twoResistanceHeatTransfer.H"
#include "addToRunTimeSelectionTable.H"
#include "generateBlendedInterfacialModels.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(massDiffusionLimitedPhaseChange, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        massDiffusionLimitedPhaseChange,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::massDiffusionLimitedPhaseChange::readCoeffs
(
    const dictionary& dict
)
{
    const phaseInterface interface(phase1_, phase2_);

    const dictionary& interfaceCompositionDict =
        dict.subDict("interfaceComposition");

    checkInterfacialModelsDict<sidedInterfaceCompositionModel>
    (
        fluid_,
        interfaceCompositionDict
    );

    interfaceCompositionModel_.reset
    (
        sidedInterfaceCompositionModel::New
        (
            interfaceCompositionDict,
            interface
        ).ptr()
    );

    const dictionary& diffusiveMassTransferDict =
        dict.subDict("diffusiveMassTransfer");

    checkBlendedInterfacialModelsDict<blendedSidedDiffusiveMassTransferModel>
    (
        fluid_,
        diffusiveMassTransferDict
    );

    diffusiveMassTransferModel_.reset
    (
        blendedSidedDiffusiveMassTransferModel::New
        (
            diffusiveMassTransferDict,
            interface,
            blendingDict<blendedSidedDiffusiveMassTransferModel>
            (
                fluid_,
                diffusiveMassTransferDict
            )
        ).ptr()
    );

    nIter_ = dict.lookupOrDefault<label>("nIter", 1);
}


Foam::wordList Foam::fv::massDiffusionLimitedPhaseChange::getSpecies() const
{
    wordHashSet species;

    forAll(phaseNames(), i)
    {
        const phaseModel& phase = i ? phase2_ : phase1_;

        if (!interfaceCompositionModel_->haveModelInThe(phase)) continue;

        species.insert(interfaceCompositionModel_->modelInThe(phase).species());
    }

    return species.toc();
}


void Foam::fv::massDiffusionLimitedPhaseChange::correctMDot() const
{
    mDot_ = dimensionedScalar(dimDensity/dimTime, 0);

    forAll(phaseNames(), i)
    {
        const phaseModel& phase = i ? phase2_ : phase1_;
        const label s = i ? +1 : -1;

        if (!interfaceCompositionModel_->haveModelInThe(phase)) continue;

        const interfaceCompositionModel& phaseIcm =
            interfaceCompositionModel_->modelInThe(phase);

        forAll(phaseIcm.species(), phaseIcmSpeciei)
        {
            const word& specieName = phaseIcm.species()[phaseIcmSpeciei];
            const label speciei = phaseIcm.thermo().species()[specieName];

            mDot_ +=
                s
               *(
                   mDotSus_[i][phaseIcmSpeciei]
                 + mDotSps_[i][phaseIcmSpeciei]*phaseIcm.thermo().Y()[speciei]
                );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::massDiffusionLimitedPhaseChange::massDiffusionLimitedPhaseChange
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phaseChange(name, modelType, mesh, dict, NullObjectRef<wordList>()),
    solver_(mesh().lookupObject<solvers::multiphaseEuler>(solver::typeName)),
    fluid_(solver_.fluid),
    phase1_(fluid_.phases()[phaseNames().first()]),
    phase2_(fluid_.phases()[phaseNames().second()]),
    interfaceCompositionModel_(nullptr),
    diffusiveMassTransferModel_(nullptr),
    nIter_(1),
    Ts_
    (
        IOobject
        (
            name + ":Ts",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        (phase1_.thermo().T()() + phase2_.thermo().T()())/2
    ),
    mDotSus_(),
    mDotSps_(),
    pressureEquationIndex_(-1),
    mDot_
    (
        IOobject
        (
            name + ":mDot",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimDensity/dimTime, 0)
    )
{
    readCoeffs(coeffs(dict));
    setSpecies(name, modelType, getSpecies());

    // Allocate the coefficient fields
    forAll(phaseNames(), i)
    {
        const phaseModel& phase = i ? phase2_ : phase1_;

        if (!interfaceCompositionModel_->haveModelInThe(phase)) continue;

        const interfaceCompositionModel& phaseIcm =
            interfaceCompositionModel_->modelInThe(phase);

        mDotSus_[i].resize(phaseIcm.species().size());

        mDotSps_[i].resize(phaseIcm.species().size());

        forAll(phaseIcm.species(), phaseIcmSpeciei)
        {
            const word& specieName = phaseIcm.species()[phaseIcmSpeciei];

            mDotSus_[i].set
            (
                phaseIcmSpeciei,
                new volScalarField::Internal
                (
                    IOobject
                    (
                        name + ":mDot" + specieName.capitalise() + "Su",
                        mesh.time().name(),
                        mesh
                    ),
                    mesh,
                    dimensionedScalar(dimDensity/dimTime, 0)
                )
            );

            mDotSps_[i].set
            (
                phaseIcmSpeciei,
                new volScalarField::Internal
                (
                    IOobject
                    (
                        name + ":mDot" + specieName.capitalise() + "Sp",
                        mesh.time().name(),
                        mesh
                    ),
                    mesh,
                    dimensionedScalar(dimDensity/dimTime, 0)
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::massDiffusionLimitedPhaseChange::Tchange() const
{
    return Ts_;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::massDiffusionLimitedPhaseChange::Lfraction() const
{
    // Heat transfer coefficients
    const Pair<tmp<volScalarField>> Hs =
        solver_.heatTransfer.Hs(phase1_, phase2_);
    const volScalarField::Internal& H1 = Hs.first();
    const volScalarField::Internal& H2 = Hs.second();

    return H2/(H1 + H2);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::massDiffusionLimitedPhaseChange::mDot() const
{
    return mDot_;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::massDiffusionLimitedPhaseChange::mDot(const label mDoti) const
{
    const word& specieName = species()[mDoti];

    tmp<volScalarField::Internal> tResult =
        volScalarField::Internal::New
        (
            name() + ":mDot" + specieName.capitalise(),
            this->mesh(),
            dimensionedScalar(dimDensity/dimTime, 0)
        );

    forAll(phaseNames(), i)
    {
        const phaseModel& phase = i ? phase2_ : phase1_;
        const label s = i ? +1 : -1;

        if (!interfaceCompositionModel_->haveModelInThe(phase)) continue;

        const interfaceCompositionModel& phaseIcm =
            interfaceCompositionModel_->modelInThe(phase);

        if (!phaseIcm.species().found(specieName)) continue;

        const label phaseIcmSpeciei = phaseIcm.species()[specieName];
        const label speciei = phaseIcm.thermo().species()[specieName];

        tResult.ref() +=
            s
           *(
               mDotSus_[i][phaseIcmSpeciei]
             + mDotSps_[i][phaseIcmSpeciei]*phaseIcm.thermo().Y()[speciei]
            );
    }

    return tResult;
}


void Foam::fv::massDiffusionLimitedPhaseChange::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    // Pressure equation (i.e., continuity, linearised in the pressure)
    if
    (
        (&alpha == &phase1_ || &alpha == &phase2_)
     && (&rho == &phase1_.rho() || &rho == &phase2_.rho())
     && &eqn.psi() == &solver_.p_rgh
    )
    {
        // Ensure that the source is up-to date if this is the first call in
        // this outer corrector
        if (pressureEquationIndex_ == 0) correctMDot();
        pressureEquationIndex_ ++;
    }

    massTransfer::addSup(alpha, rho, eqn);
}


void Foam::fv::massDiffusionLimitedPhaseChange::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& heOrYi,
    fvMatrix<scalar>& eqn
) const
{
    const label index = this->index(phaseNames(), alpha.group());

    const ThermoRefPair<multicomponentThermo>& mcThermos =
        thermos().thermos<multicomponentThermo>();

    const word specieName = heOrYi.member();

    // Mass fraction equation
    if (mcThermos.valid()[index] && mcThermos[index].containsSpecie(specieName))
    {
        // A non-transferring specie. Do not add a source.
        if (!species().found(specieName)) return;

        // A transferring specie. Add a (potentially) linearised source.
        forAll(phaseNames(), i)
        {
            const phaseModel& phase = i ? phase2_ : phase1_;

            if (!interfaceCompositionModel_->haveModelInThe(phase)) continue;

            const interfaceCompositionModel& phaseIcm =
                interfaceCompositionModel_->modelInThe(phase);

            if (!phaseIcm.species().found(specieName)) continue;

            const label phaseIcmSpeciei = phaseIcm.species()[specieName];
            const label speciei = phaseIcm.thermo().species()[specieName];

            if (i == index)
            {
                eqn +=
                    mDotSus_[i][phaseIcmSpeciei]
                  + fvm::Sp(mDotSps_[i][phaseIcmSpeciei], eqn.psi());
            }
            else
            {
                eqn -=
                    mDotSus_[i][phaseIcmSpeciei]
                  + mDotSps_[i][phaseIcmSpeciei]*phaseIcm.thermo().Y()[speciei];
            }
        }

        return;
    }

    // Something else. Fall back.
    phaseChange::addSup(alpha, rho, heOrYi, eqn);
}


void Foam::fv::massDiffusionLimitedPhaseChange::correct()
{
    Info<< type() << ": " << name() << endl << incrIndent;

    const volScalarField::Internal& T1 = phase1_.thermo().T();
    const volScalarField::Internal& T2 = phase2_.thermo().T();

    // Heat transfer coefficients
    const Pair<tmp<volScalarField>> Hs =
        solver_.heatTransfer.Hs(phase1_, phase2_, scalar(0));
    const volScalarField::Internal& H1 = Hs.first();
    const volScalarField::Internal& H2 = Hs.second();

    // Stabilisation heat transfer coefficient
    static const dimensionedScalar HSmall
    (
        "small",
        heatTransferModel::dimK,
        small
    );

    // Mass transfer coefficients
    Pair<autoPtr<volScalarField::Internal>> KPtrs;
    forAll(phaseNames(), i)
    {
        const phaseModel& phase = i ? phase2_ : phase1_;

        if (!interfaceCompositionModel_->haveModelInThe(phase)) continue;

        KPtrs[i].set(diffusiveMassTransferModel_->KinThe(phase).ptr());
    }

    // Mass diffusivities
    Pair<PtrList<volScalarField::Internal>> Ds;
    forAll(phaseNames(), i)
    {
        const phaseModel& phase = i ? phase2_ : phase1_;

        if (!interfaceCompositionModel_->haveModelInThe(phase)) continue;

        interfaceCompositionModel& phaseIcm =
            interfaceCompositionModel_->modelInThe(phase);

        Ds[i].setSize(phaseIcm.species().size());

        forAll(phaseIcm.species(), phaseIcmSpeciei)
        {
            const word& specieName = phaseIcm.species()[phaseIcmSpeciei];

            Ds[i].set
            (
                phaseIcmSpeciei,
                phaseIcm.D(specieName).ptr()
            );
        }
    }

    // Iterative solution of the interface heat-mass transfer balance
    for (label iteri = 0; iteri < nIter_; ++ iteri)
    {
        tmp<volScalarField::Internal> mDotL =
            volScalarField::Internal::New
            (
                name() + ":mDotL",
                this->mesh(),
                dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
            );
        tmp<volScalarField::Internal> mDotLPrime =
            volScalarField::Internal::New
            (
                name() + ":mDotLPrime",
                this->mesh(),
                dimensionedScalar(mDotL().dimensions()/dimTemperature, 0)
            );

        // Add latent heats from forward and backward models
        forAll(phaseNames(), i)
        {
            const phaseModel& phase = i ? phase2_ : phase1_;
            const label s = i ? +1 : -1;

            if (!interfaceCompositionModel_->haveModelInThe(phase)) continue;

            interfaceCompositionModel& phaseIcm =
                interfaceCompositionModel_->modelInThe(phase);

            forAll(phaseIcm.species(), phaseIcmSpeciei)
            {
                const word& specieName = phaseIcm.species()[phaseIcmSpeciei];
                const label speciei = phaseIcm.thermo().species()[specieName];
                const label mDoti = species()[specieName];

                const volScalarField::Internal Yf
                (
                    phaseIcm.Yf(specieName, vifToVf(Ts_))
                );
                const volScalarField::Internal YfPrime
                (
                    phaseIcm.YfPrime(specieName, vifToVf(Ts_))
                );

                const volScalarField::Internal rhoKDL
                (
                    phase.rho()()
                   *KPtrs[i]()
                   *Ds[i][phaseIcmSpeciei]
                   *s*L(Ts_, mDoti)
                );

                mDotL.ref() += rhoKDL*(Yf - phaseIcm.thermo().Y()[speciei]);
                mDotLPrime.ref() += rhoKDL*YfPrime;
            }
        }

        // Update the interface temperature by applying one step of newton's
        // method to the interface relation
        Ts_ -=
            (H1*(Ts_ - T1) + H2*(Ts_ - T2) + mDotL)
           /(max(H1 + H2 + mDotLPrime, HSmall));

        Info<< indent << "min/mean/max Ts"
            << " = " << gMin(Ts_.primitiveField())
            << '/' << gAverage(Ts_.primitiveField())
            << '/' << gMax(Ts_.primitiveField())
            << endl;

        // Update the interface compositions
        forAll(phaseNames(), i)
        {
            const phaseModel& phase = i ? phase2_ : phase1_;

            if (!interfaceCompositionModel_->haveModelInThe(phase)) continue;

            interfaceCompositionModel& phaseIcm =
                interfaceCompositionModel_->modelInThe(phase);

            phaseIcm.update(vifToVf(Ts_));
        }
    }

    // Build the coefficients of the phase change rate
    forAll(phaseNames(), i)
    {
        const phaseModel& phase = i ? phase2_ : phase1_;

        if (!interfaceCompositionModel_->haveModelInThe(phase)) continue;

        interfaceCompositionModel& phaseIcm =
            interfaceCompositionModel_->modelInThe(phase);

        forAll(phaseIcm.species(), phaseIcmSpeciei)
        {
            const word& specieName = phaseIcm.species()[phaseIcmSpeciei];

            const volScalarField::Internal KD
            (
                KPtrs[i]()*Ds[i][phaseIcmSpeciei]
            );

            const volScalarField::Internal Yf
            (
                phaseIcm.Yf(specieName, vifToVf(Ts_))
            );

            mDotSus_[i][phaseIcmSpeciei] = phase.rho()()*KD*Yf;
            mDotSps_[i][phaseIcmSpeciei] = -phase.rho()()*KD;
        }
    }

    // Reset the p_rgh equation solution counter
    pressureEquationIndex_ = 0;

    // Correct the total phase change rate
    correctMDot();

    Info<< indent << "min/mean/max mDot"
        << " = " << gMin(mDot_.primitiveField())
        << '/' << gAverage(mDot_.primitiveField())
        << '/' << gMax(mDot_.primitiveField())
        << endl;

    Info<< decrIndent;
}


bool Foam::fv::massDiffusionLimitedPhaseChange::read(const dictionary& dict)
{
    if (phaseChange::read(dict))
    {
        readCoeffs(coeffs(dict));
        reSetSpecies(getSpecies());
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
