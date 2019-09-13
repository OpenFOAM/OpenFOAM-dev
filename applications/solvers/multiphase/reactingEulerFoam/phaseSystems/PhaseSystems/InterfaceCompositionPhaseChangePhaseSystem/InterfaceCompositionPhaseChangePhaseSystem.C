/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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

#include "InterfaceCompositionPhaseChangePhaseSystem.H"
#include "interfaceCompositionModel.H"
#include "heatTransferModel.H"
#include "diffusiveMassTransferModel.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::dmdtTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
iDmdts() const
{
    autoPtr<phaseSystem::dmdtTable> iDmdtsPtr(new phaseSystem::dmdtTable);

    phaseSystem::dmdtTable& iDmdts = iDmdtsPtr();

    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];

        forAllConstIter(phasePair, pair, pairIter)
        {
            const autoPtr<interfaceCompositionModel>& compositionModelPtr =
                interfaceCompositionModelIter()[pairIter.index()];

            if (!compositionModelPtr.valid()) continue;

            const interfaceCompositionModel& compositionModel =
                compositionModelPtr();

            const phaseModel& phase = *pairIter;

            forAllConstIter
            (
                hashedWordList,
                compositionModel.species(),
                memberIter
            )
            {
                const word& member = *memberIter;

                const label dmidtSign = pairIter.index() == 0 ? +1 : -1;

                tmp<volScalarField> dmidt
                (
                    *(*iDmdtSu_[pair])[member]
                  + *(*iDmdtSp_[pair])[member]*phase.Y(member)
                );

                if (iDmdts.found(pair))
                {
                    *iDmdts[pair] += dmidtSign*dmidt;
                }
                else
                {
                    iDmdts.insert(pair, (dmidtSign*dmidt).ptr());
                }
            }
        }
    }

    return iDmdtsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::dmidtTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
iDmidts() const
{
    autoPtr<phaseSystem::dmidtTable> iDmidtsPtr(new phaseSystem::dmidtTable);

    phaseSystem::dmidtTable& iDmidts = iDmidtsPtr();

    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];

        if (!iDmidts.found(pair))
        {
            iDmidts.insert(pair, new HashPtrTable<volScalarField>());
        }

        forAllConstIter(phasePair, pair, pairIter)
        {
            const autoPtr<interfaceCompositionModel>& compositionModelPtr =
                interfaceCompositionModelIter()[pairIter.index()];

            if (!compositionModelPtr.valid()) continue;

            const interfaceCompositionModel& compositionModel =
                compositionModelPtr();

            const phaseModel& phase = *pairIter;

            forAllConstIter
            (
                hashedWordList,
                compositionModel.species(),
                memberIter
            )
            {
                const word& member = *memberIter;

                const label dmidtSign = pairIter.index() == 0 ? +1 : -1;

                tmp<volScalarField> dmidt
                (
                    *(*iDmdtSu_[pair])[member]
                  + *(*iDmdtSp_[pair])[member]*phase.Y(member)
                );

                if (iDmidts[pair]->found(member))
                {
                    *(*iDmidts[pair])[member] += dmidtSign*dmidt;
                }
                else
                {
                    iDmidts[pair]->insert(member, (dmidtSign*dmidt).ptr());
                }
            }
        }
    }

    return iDmidtsPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
InterfaceCompositionPhaseChangePhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh),
    nInterfaceCorrectors_
    (
        this->template lookupOrDefault<label>("nInterfaceCorrectors", 1)
    )
{
    this->generatePairsAndSubModels
    (
        "interfaceComposition",
        interfaceCompositionModels_
    );

    this->generatePairsAndSubModels
    (
        "diffusiveMassTransfer",
        diffusiveMassTransferModels_,
        false
    );

    // Check that models have been specified in the correct combinations
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];

        if (!this->diffusiveMassTransferModels_.found(pair))
        {
            FatalErrorInFunction
                << "A diffusive mass transfer model the " << pair
                << " pair is not specified. This is required by the "
                << "corresponding interface composition model."
                << exit(FatalError);
        }

        forAllConstIter(phasePair, pair, pairIter)
        {
            if
            (
                interfaceCompositionModelIter()[pairIter.index()].valid()
             && !diffusiveMassTransferModels_[pair][pairIter.index()].valid()
            )
            {
                FatalErrorInFunction
                    << "A mass transfer model for the " << (*pairIter).name()
                    << " side of the " << pair << " pair is not "
                    << "specified. This is required by the corresponding "
                    << "interface composition model."
                    << exit(FatalError);
            }
        }

        if
        (
            !this->heatTransferModels_.found(pair)
         || !this->heatTransferModels_[pair].first().valid()
         || !this->heatTransferModels_[pair].second().valid()
        )
        {
             FatalErrorInFunction
                 << "A heat transfer model for both sides of the " << pair
                 << "pair is not specified. This is required by the "
                 << "corresponding interface composition model"
                 << exit(FatalError);
        }
    }

    // Generate mass transfer fields, initially assumed to be zero
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];

        iDmdtSu_.insert(pair, new HashPtrTable<volScalarField>());
        iDmdtSp_.insert(pair, new HashPtrTable<volScalarField>());

        forAllConstIter(phasePair, pair, pairIter)
        {
            const autoPtr<interfaceCompositionModel>& compositionModelPtr =
                interfaceCompositionModelIter()[pairIter.index()];

            if (!compositionModelPtr.valid()) continue;

            const interfaceCompositionModel& compositionModel =
                compositionModelPtr();

            forAllConstIter
            (
                hashedWordList,
                compositionModel.species(),
                memberIter
            )
            {
                const word& member = *memberIter;

                iDmdtSu_[pair]->insert
                (
                    member,
                    zeroVolField<scalar>(pair, "iDmdtSu", dimDensity/dimTime)
                   .ptr()
                );

                iDmdtSp_[pair]->insert
                (
                    member,
                    zeroVolField<scalar>(pair, "iDmdtSp", dimDensity/dimTime)
                   .ptr()
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
~InterfaceCompositionPhaseChangePhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    autoPtr<phaseSystem::dmdtTable> iDmdtsPtr = this->iDmdts();

    const phaseSystem::dmdtTable& iDmdts = iDmdtsPtr();

    forAllConstIter(phaseSystem::dmdtTable, iDmdts, iDmdtIter)
    {
        const phasePair& pair = this->phasePairs_[iDmdtIter.key()];
        const phaseModel& phase = pair.phase1();
        const phaseModel& otherPhase = pair.phase2();

        addField(phase, "dmdt", *iDmdtIter(), dmdts);
        addField(otherPhase, "dmdt", - *iDmdtIter(), dmdts);
    }

    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
momentumTransfer()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransfer();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtU(iDmdts(), eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
momentumTransferf()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransferf();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtU(iDmdts(), eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    this->addDmidtHe(iDmidts(), eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::specieTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
specieTransfer() const
{
    autoPtr<phaseSystem::specieTransferTable> eqnsPtr =
        BasePhaseSystem::specieTransfer();

    phaseSystem::specieTransferTable& eqns = eqnsPtr();

    // Explicit
    /*
    this->addDmidtY(iDmidts(), eqns);
    */

    // Semi-implicit
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];

        forAllConstIter(phasePair, pair, pairIter)
        {
            const autoPtr<interfaceCompositionModel>& compositionModelPtr =
                interfaceCompositionModelIter()[pairIter.index()];

            if (!compositionModelPtr.valid()) continue;

            const interfaceCompositionModel& compositionModel =
                compositionModelPtr();

            const phaseModel& phase = *pairIter;
            const phaseModel& otherPhase = pairIter.otherPhase();

            forAllConstIter
            (
                hashedWordList,
                compositionModel.species(),
                memberIter
            )
            {
                const word& member = *memberIter;

                // Implicit transport through this phase
                *eqns[phase.Y(member).name()] +=
                    *(*iDmdtSu_[pair])[member]
                  + fvm::Sp(*(*iDmdtSp_[pair])[member], phase.Y(member));

                // Explicit transport out of the other phase
                if (eqns.found(IOobject::groupName(member, otherPhase.name())))
                {
                    *eqns[otherPhase.Y(member).name()] -=
                        *(*iDmdtSu_[pair])[member]
                      + *(*iDmdtSp_[pair])[member]*phase.Y(member);
                }
            }
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
correct()
{
    BasePhaseSystem::correct();

    // Sum up the contribution from each interface composition model
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];

        const volScalarField& Tf(*this->Tf_[pair]);

        forAllConstIter(phasePair, pair, pairIter)
        {
            const autoPtr<interfaceCompositionModel>& compositionModelPtr =
                interfaceCompositionModelIter()[pairIter.index()];

            if (!compositionModelPtr.valid()) continue;

            const interfaceCompositionModel& compositionModel =
                compositionModelPtr();

            const phaseModel& phase = *pairIter;

            const volScalarField K
            (
                diffusiveMassTransferModels_[pair][pairIter.index()]->K()
            );

            forAllConstIter
            (
                hashedWordList,
                compositionModel.species(),
                memberIter
            )
            {
                const word& member = *memberIter;

                const volScalarField KD(K*compositionModel.D(member));
                const volScalarField Yf(compositionModel.Yf(member, Tf));

                *(*iDmdtSu_[pair])[member] = phase.rho()*KD*Yf;
                *(*iDmdtSp_[pair])[member] = - phase.rho()*KD;
            }
        }
    }
}


template<class BasePhaseSystem>
void Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
correctInterfaceThermo()
{
    // This loop solves for the interface temperatures, Tf, and updates the
    // interface composition models.
    //
    // In the presence of thermally-coupled mass transfer, the relation between
    // heat transfers across the interface between phases 1 and 2 is:
    //
    //                         Q1 + Q2 = mDot*L
    //     H1*(Tf - T1) + H2*(Tf - T1) = K*rho*(Yfi - Yi)*Li
    //
    // Where Q1 and Q2 are the net transfer into phases 1 and 2 respectively,
    // H1 and H2 are the heat transfer coefficients on either side, Tf is the
    // temperature at the interface, mDot is the mass transfer rate from phase
    // 2 to phase 1, and L is the latent heat of phase 2 minus phase 1, K is
    // the diffusive mass transfer coefficient, Yfi - Yi is the concentration
    // difference of a transferring specie between the interface and the bulk
    // driving the transfer, Li is the latent heat change of the specie, and
    // rho is the density in the phase in which the diffusive mass transfer is
    // being represented.
    //
    // Yfi is likely to be a strong non-linear (typically exponential) function
    // of Tf, so the solution for the temperature is newton-accelerated.

    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];

        const volScalarField H1(this->heatTransferModels_[pair].first()->K());
        const volScalarField H2(this->heatTransferModels_[pair].second()->K());
        const dimensionedScalar HSmall("small", heatTransferModel::dimK, small);

        volScalarField& Tf = *this->Tf_[pair];

        for (label i = 0; i < nInterfaceCorrectors_; ++ i)
        {
            tmp<volScalarField> mDotL =
                zeroVolField<scalar>
                (
                    pair,
                    "mDotL",
                    dimEnergy/dimVolume/dimTime
                );
            tmp<volScalarField> mDotLPrime =
                zeroVolField<scalar>
                (
                    pair,
                    "mDotLPrime",
                    mDotL().dimensions()/dimTemperature
                );

            // Add latent heats from forward and backward models
            if (this->interfaceCompositionModels_[pair].first().valid())
            {
                this->interfaceCompositionModels_[pair].first()->addMDotL
                (
                    diffusiveMassTransferModels_[pair].first()->K(),
                    Tf,
                    mDotL.ref(),
                    mDotLPrime.ref()
                );
            }
            if (this->interfaceCompositionModels_[pair].second().valid())
            {
                this->interfaceCompositionModels_[pair].second()->addMDotL
                (
                  - diffusiveMassTransferModels_[pair].second()->K(),
                    Tf,
                    mDotL.ref(),
                    mDotLPrime.ref()
                );
            }

            // Update the interface temperature by applying one step of newton's
            // method to the interface relation
            Tf -=
                (
                    H1*(Tf - pair.phase1().thermo().T())
                  + H2*(Tf - pair.phase2().thermo().T())
                  - mDotL
                )
               /(
                    max(H1 + H2 - mDotLPrime, HSmall)
                );

            Tf.correctBoundaryConditions();

            Info<< "Tf." << pair.name()
                << ": min = " << min(Tf.primitiveField())
                << ", mean = " << average(Tf.primitiveField())
                << ", max = " << max(Tf.primitiveField())
                << endl;

            // Update the interface compositions
            if (this->interfaceCompositionModels_[pair].first().valid())
            {
                this->interfaceCompositionModels_[pair].first()->update(Tf);
            }
            if (this->interfaceCompositionModels_[pair].second().valid())
            {
                this->interfaceCompositionModels_[pair].second()->update(Tf);
            }
        }
    }
}


template<class BasePhaseSystem>
bool Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
