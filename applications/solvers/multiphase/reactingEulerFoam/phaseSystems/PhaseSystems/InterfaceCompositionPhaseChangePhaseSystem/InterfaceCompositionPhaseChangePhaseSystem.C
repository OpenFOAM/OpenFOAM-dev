/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
InterfaceCompositionPhaseChangePhaseSystem
(
    const fvMesh& mesh
)
:
    HeatAndMassTransferPhaseSystem<BasePhaseSystem>(mesh)
{
    this->generatePairsAndSubModels
    (
        "interfaceComposition",
        interfaceCompositionModels_
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
~InterfaceCompositionPhaseChangePhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::massTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
massTransfer() const
{
    // Create a mass transfer matrix for each species of each phase
    autoPtr<phaseSystem::massTransferTable> eqnsPtr
    (
        new phaseSystem::massTransferTable()
    );

    phaseSystem::massTransferTable& eqns = eqnsPtr();

    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        const PtrList<volScalarField>& Yi = phase.Y();

        forAll(Yi, i)
        {
            eqns.insert
            (
                Yi[i].name(),
                new fvScalarMatrix(Yi[i], dimMass/dimTime)
            );
        }
    }

    // Reset the interfacial mass flow rates
    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const phasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        *this->dmdt_[pair] =
            *this->dmdtExplicit_[pair];

        *this->dmdtExplicit_[pair] =
            dimensionedScalar("zero", dimDensity/dimTime, 0);
    }

    // Sum up the contribution from each interface composition model
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const interfaceCompositionModel& compositionModel
        (
            interfaceCompositionModelIter()
        );

        const phasePair& pair
        (
            this->phasePairs_[interfaceCompositionModelIter.key()]
        );
        const phaseModel& phase = pair.phase1();
        const phaseModel& otherPhase = pair.phase2();
        const phasePairKey key(phase.name(), otherPhase.name());

        const volScalarField& Tf(*this->Tf_[key]);

        volScalarField& dmdtExplicit(*this->dmdtExplicit_[key]);
        volScalarField& dmdt(*this->dmdt_[key]);

        scalar dmdtSign(Pair<word>::compare(this->dmdt_.find(key).key(), key));

        const volScalarField K
        (
            this->massTransferModels_[key][phase.name()]->K()
        );

        forAllConstIter
        (
            hashedWordList,
            compositionModel.species(),
            memberIter
        )
        {
            const word& member = *memberIter;

            const word name
            (
                IOobject::groupName(member, phase.name())
            );

            const word otherName
            (
                IOobject::groupName(member, otherPhase.name())
            );

            const volScalarField KD
            (
                K*compositionModel.D(member)
            );

            const volScalarField Yf
            (
                compositionModel.Yf(member, Tf)
            );

            // Implicit transport through the phase
            *eqns[name] +=
                phase.rho()*KD*Yf
              - fvm::Sp(phase.rho()*KD, eqns[name]->psi());

            // Sum the mass transfer rate
            dmdtExplicit += dmdtSign*phase.rho()*KD*Yf;
            dmdt -= dmdtSign*phase.rho()*KD*eqns[name]->psi();

            // Explicit transport out of the other phase
            if (eqns.found(otherName))
            {
                *eqns[otherName] -=
                    otherPhase.rho()*KD*compositionModel.dY(member, Tf);
            }
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
correctThermo()
{
    BasePhaseSystem::correctThermo();

    // This loop solves for the interface temperatures, Tf, and updates the
    // interface composition models.
    //
    // The rate of heat transfer to the interface must equal the latent heat
    // consumed at the interface, i.e.:
    //
    // H1*(T1 - Tf) + H2*(T2 - Tf) == mDotL
    //                             == K*rho*(Yfi - Yi)*Li
    //
    // Yfi is likely to be a strong non-linear (typically exponential) function
    // of Tf, so the solution for the temperature is newton-accelerated

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const phasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        const phasePairKey key12(pair.first(), pair.second(), true);
        const phasePairKey key21(pair.second(), pair.first(), true);

        volScalarField H1(this->heatTransferModels_[pair][pair.first()]->K());
        volScalarField H2(this->heatTransferModels_[pair][pair.second()]->K());
        dimensionedScalar HSmall("small", heatTransferModel::dimK, SMALL);

        volScalarField mDotL
        (
            IOobject
            (
                "mDotL",
                this->mesh().time().timeName(),
                this->mesh()
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0)
        );
        volScalarField mDotLPrime
        (
            IOobject
            (
                "mDotLPrime",
                this->mesh().time().timeName(),
                this->mesh()
            ),
            this->mesh(),
            dimensionedScalar("zero", mDotL.dimensions()/dimTemperature, 0)
        );

        volScalarField& Tf = *this->Tf_[pair];

        // Add latent heats from forward and backward models
        if (this->interfaceCompositionModels_.found(key12))
        {
            this->interfaceCompositionModels_[key12]->addMDotL
            (
                this->massTransferModels_[pair][pair.first()]->K(),
                Tf,
                mDotL,
                mDotLPrime
            );
        }
        if (this->interfaceCompositionModels_.found(key21))
        {
            this->interfaceCompositionModels_[key21]->addMDotL
            (
                this->massTransferModels_[pair][pair.second()]->K(),
                Tf,
                mDotL,
                mDotLPrime
            );
        }

        // Update the interface temperature by applying one step of newton's
        // method to the interface relation
        Tf -=
            (
                H1*(Tf - pair.phase1().thermo().T())
              + H2*(Tf - pair.phase2().thermo().T())
              + mDotL
            )
           /(
                max(H1 + H2 + mDotLPrime, HSmall)
            );

        Tf.correctBoundaryConditions();

        Info<< "Tf." << pair.name()
            << ": min = " << min(Tf.internalField())
            << ", mean = " << average(Tf.internalField())
            << ", max = " << max(Tf.internalField())
            << endl;

        // Update the interface compositions
        if (this->interfaceCompositionModels_.found(key12))
        {
            this->interfaceCompositionModels_[key12]->update(Tf);
        }
        if (this->interfaceCompositionModels_.found(key21))
        {
            this->interfaceCompositionModels_[key21]->update(Tf);
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
