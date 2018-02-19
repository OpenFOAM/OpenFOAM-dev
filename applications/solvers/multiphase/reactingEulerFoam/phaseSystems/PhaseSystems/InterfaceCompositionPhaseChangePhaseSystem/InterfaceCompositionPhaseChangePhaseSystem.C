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

#include "InterfaceCompositionPhaseChangePhaseSystem.H"
#include "interfaceCompositionModel.H"


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
addMomentumTransfer(phaseSystem::momentumTransferTable& eqns) const
{
    // Source term due to mass transfer
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

        const volVectorField& U1(pair.phase1().U());
        const volVectorField& U2(pair.phase2().U());

        const volScalarField dmdt(this->iDmdt(pair));
        const volScalarField dmdt21(posPart(dmdt));
        const volScalarField dmdt12(negPart(dmdt));

        *eqns[pair.phase1().name()] += dmdt21*U2 - fvm::Sp(dmdt21, U1);
        *eqns[pair.phase2().name()] -= dmdt12*U1 - fvm::Sp(dmdt12, U2);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
InterfaceCompositionPhaseChangePhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generatePairsAndSubModels
    (
        "interfaceComposition",
        interfaceCompositionModels_
    );

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

        // Initially assume no mass transfer
        iDmdt_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("iDmdt", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimDensity/dimTime, 0)
            )
        );

        iDmdtExplicit_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("iDmdtExplicit", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar("zero", dimDensity/dimTime, 0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
~InterfaceCompositionPhaseChangePhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::iDmdt
(
    const phasePairKey& key
) const
{
    const scalar dmdtSign(Pair<word>::compare(iDmdt_.find(key).key(), key));

    return dmdtSign**iDmdt_[key];
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::iDmdt
(
    const Foam::phaseModel& phase
) const
{
    tmp<volScalarField> tiDmdt
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("iDmdt", phase.name()),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("zero", dimDensity/dimTime, 0)
        )
    );

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

        if (pair.contains(phase))
        {
            tiDmdt.ref() += this->iDmdt
            (
                phasePairKey(phase.name(), pair.otherPhase(phase).name(), false)
            );
        }
    }

    return tiDmdt;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::dmdt
(
    const Foam::phaseModel& phase
) const
{
    tmp<volScalarField> tDmdt
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("dmdt", phase.name()),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar("zero", dimDensity/dimTime, 0)
        )
    );

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

        if (pair.contains(phase))
        {
            tDmdt.ref() += this->dmdt
            (
                phasePairKey(phase.name(), pair.otherPhase(phase).name(), false)
            );
        }
    }

    return tDmdt;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
momentumTransfer() const
{
    autoPtr<phaseSystem::momentumTransferTable>
        eqnsPtr(BasePhaseSystem::momentumTransfer());

    addMomentumTransfer(eqnsPtr());

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
momentumTransferf() const
{
    autoPtr<phaseSystem::momentumTransferTable>
        eqnsPtr(BasePhaseSystem::momentumTransferf());

    addMomentumTransfer(eqnsPtr());

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

    // Source term due to mass transfer
    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        if
        (
            this->heatTransferModels_.found(phasePairIter.key())
        )
        {
            const phasePair& pair(phasePairIter());

            if (pair.ordered())
            {
                continue;
            }

            const phaseModel& phase1 = pair.phase1();
            const phaseModel& phase2 = pair.phase2();

            const volScalarField& he1(phase1.thermo().he());
            const volScalarField& he2(phase2.thermo().he());

            const volScalarField& K1(phase1.K());
            const volScalarField& K2(phase2.K());

            const volScalarField dmdt(this->dmdt(pair));
            const volScalarField dmdt21(posPart(dmdt));
            const volScalarField dmdt12(negPart(dmdt));
            const volScalarField& Tf(*this->Tf_[pair]);

            *eqns[phase1.name()] +=
                dmdt21*(phase1.thermo().he(phase1.thermo().p(), Tf))
              - fvm::Sp(dmdt21, he1)
              + dmdt21*(K2 - K1);

            *eqns[phase2.name()] -=
                dmdt12*(phase2.thermo().he(phase2.thermo().p(), Tf))
              - fvm::Sp(dmdt12, he2)
              + dmdt12*(K1 - K2);
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::massTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
massTransfer() const
{
    autoPtr<phaseSystem::massTransferTable> eqnsPtr =
        BasePhaseSystem::massTransfer();

    phaseSystem::massTransferTable& eqns = eqnsPtr();

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

        *this->iDmdt_[pair] =
            *this->iDmdtExplicit_[pair];

        *this->iDmdtExplicit_[pair] =
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

        volScalarField& iDmdtExplicit_(*this->iDmdtExplicit_[key]);
        volScalarField& iDmdt_(*this->iDmdt_[key]);

        scalar dmdtSign(Pair<word>::compare(this->iDmdt_.find(key).key(), key));

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
            iDmdtExplicit_ += dmdtSign*phase.rho()*KD*Yf;
            iDmdt_ -= dmdtSign*phase.rho()*KD*eqns[name]->psi();

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
    phaseSystem::correctThermo();

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
        dimensionedScalar HSmall("small", heatTransferModel::dimK, small);

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
            << ": min = " << min(Tf.primitiveField())
            << ", mean = " << average(Tf.primitiveField())
            << ", max = " << max(Tf.primitiveField())
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
