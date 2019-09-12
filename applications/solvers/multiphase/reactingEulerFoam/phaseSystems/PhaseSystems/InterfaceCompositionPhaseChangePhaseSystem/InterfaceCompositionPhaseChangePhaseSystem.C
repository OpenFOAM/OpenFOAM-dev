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
#include "massTransferModel.H"
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
        const interfaceCompositionModel& compositionModel =
            interfaceCompositionModelIter();

        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];
        const phaseModel& phase = pair.phase1();
        const phaseModel& otherPhase = pair.phase2();
        const phasePair& unorderedPair =
            this->phasePairs_[phasePair(phase, otherPhase)];

        forAllConstIter(hashedWordList, compositionModel.species(), memberIter)
        {
            const word& member = *memberIter;

            tmp<volScalarField> dmidt
            (
                Pair<word>::compare(pair, unorderedPair)
               *(
                    *(*iDmdtSu_[pair])[member]
                  + *(*iDmdtSp_[pair])[member]*phase.Y(member)
                )
            );

            if (iDmdts.found(unorderedPair))
            {
                *iDmdts[unorderedPair] += dmidt;
            }
            else
            {
                iDmdts.insert(unorderedPair, dmidt.ptr());
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
        const interfaceCompositionModel& compositionModel =
            interfaceCompositionModelIter();

        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];
        const phaseModel& phase = pair.phase1();
        const phaseModel& otherPhase = pair.phase2();
        const phasePair& unorderedPair =
            this->phasePairs_[phasePair(phase, otherPhase)];

        if (!iDmidts.found(unorderedPair))
        {
            iDmidts.insert(unorderedPair, new HashPtrTable<volScalarField>());
        }

        forAllConstIter(hashedWordList, compositionModel.species(), memberIter)
        {
            const word& member = *memberIter;

            tmp<volScalarField> dmidt
            (
                Pair<word>::compare(pair, unorderedPair)
               *(
                    *(*iDmdtSu_[pair])[member]
                  + *(*iDmdtSp_[pair])[member]*phase.Y(member)
                )
            );

            if (iDmidts[unorderedPair]->found(member))
            {
                *(*iDmidts[unorderedPair])[member] += dmidt;
            }
            else
            {
                iDmidts[unorderedPair]->insert(member, dmidt.ptr());
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
        "massTransfer",
        massTransferModels_,
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
        const phaseModel& phase = pair.phase1();
        const phaseModel& otherPhase = pair.phase2();

        if (!pair.ordered())
        {
            FatalErrorInFunction
                << "An interfacial composition model is specified for the "
                << "unordered " << pair << " pair. Composition models only "
                << "apply to ordered pairs. An entry for a "
                << phasePairKey("A", "B", true) << " pair means a model for "
                << "the A side of the A-B interface; i.e., \"A in the presence "
                << "of B\""
                << exit(FatalError);
        }


        const phasePairKey key(phase.name(), otherPhase.name());

        if (!this->phasePairs_.found(key))
        {
            FatalErrorInFunction
                << "A mass transfer model the " << key << " pair is not "
                << "specified. This is required by the corresponding interface "
                << "composition model."
                << exit(FatalError);
        }

        const phasePair& uoPair = this->phasePairs_[key];

        if (!massTransferModels_[uoPair][uoPair.index(phase)].valid())
        {
            FatalErrorInFunction
                << "A mass transfer model for the " << pair.phase1().name()
                << " side of the " << uoPair << " pair is not "
                << "specified. This is required by the corresponding interface "
                << "composition model."
                << exit(FatalError);
        }
    }
    forAllConstIter
    (
        massTransferModelTable,
        massTransferModels_,
        massTransferModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[massTransferModelIter.key()];

        if (!this->heatTransferModels_.found(pair))
        {
             FatalErrorInFunction
                 << "A heat transfer model for " << pair << " pair is not "
                 << "specified. This is required by the corresponding species "
                 << "transfer model"
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
        const interfaceCompositionModel& compositionModel =
            interfaceCompositionModelIter();

        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];

        iDmdtSu_.insert(pair, new HashPtrTable<volScalarField>());
        iDmdtSp_.insert(pair, new HashPtrTable<volScalarField>());

        forAllConstIter(hashedWordList, compositionModel.species(), memberIter)
        {
            const word& member = *memberIter;

            iDmdtSu_[pair]->insert
            (
                member,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("iDmdtSu", pair.name()),
                        this->mesh().time().timeName(),
                        this->mesh()
                    ),
                    this->mesh(),
                    dimensionedScalar(dimDensity/dimTime, 0)
                )
            );

            iDmdtSp_[pair]->insert
            (
                member,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("iDmdtSp", pair.name()),
                        this->mesh().time().timeName(),
                        this->mesh()
                    ),
                    this->mesh(),
                    dimensionedScalar(dimDensity/dimTime, 0)
                )
            );
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
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    NotImplemented;

    return phaseSystem::dmdt(key);
}


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

        this->addField(phase, "dmdt", *iDmdtIter(), dmdts);
        this->addField(otherPhase, "dmdt", - *iDmdtIter(), dmdts);
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
Foam::autoPtr<Foam::phaseSystem::massTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
massTransfer() const
{
    autoPtr<phaseSystem::massTransferTable> eqnsPtr =
        BasePhaseSystem::massTransfer();

    phaseSystem::massTransferTable& eqns = eqnsPtr();

    // Explicit
    /*
    this->addDmidtY(iDmidts(), eqnsPtr());
    */

    // Semi-implicit
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const interfaceCompositionModel& compositionModel =
            interfaceCompositionModelIter();

        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];
        const phaseModel& phase = pair.phase1();
        const phaseModel& otherPhase = pair.phase2();

        forAllConstIter(hashedWordList, compositionModel.species(), memberIter)
        {
            const word& member = *memberIter;
            const word name = IOobject::groupName(member, phase.name());
            const word otherName =
                IOobject::groupName(member, otherPhase.name());

            // Implicit transport through this phase
            *eqns[name] +=
                *(*iDmdtSu_[pair])[member]
              + fvm::Sp(*(*iDmdtSp_[pair])[member], phase.Y(member));

            // Explicit transport out of the other phase
            if (eqns.found(otherName))
            {
                *eqns[otherName] -=
                    *(*iDmdtSu_[pair])[member]
                  + *(*iDmdtSp_[pair])[member]*phase.Y(member);
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
        const interfaceCompositionModel& compositionModel =
            interfaceCompositionModelIter();

        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];
        const phaseModel& phase = pair.phase1();
        const phaseModel& otherPhase = pair.phase2();
        const phasePair& unorderedPair =
            this->phasePairs_[phasePair(phase, otherPhase)];

        const volScalarField& Tf(*this->Tf_[unorderedPair]);

        const volScalarField K
        (
            massTransferModels_[unorderedPair][unorderedPair.index(phase)]->K()
        );

        forAllConstIter(hashedWordList, compositionModel.species(), memberIter)
        {
            const word& member = *memberIter;

            const volScalarField KD(K*compositionModel.D(member));
            const volScalarField Yf(compositionModel.Yf(member, Tf));

            *(*iDmdtSu_[pair])[member] = phase.rho()*KD*Yf;
            *(*iDmdtSp_[pair])[member] = - phase.rho()*KD;
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
        typename BasePhaseSystem::heatTransferModelTable,
        this->heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[heatTransferModelIter.key()];

        const phasePairKey key12(pair.first(), pair.second(), true);
        const phasePairKey key21(pair.second(), pair.first(), true);

        const volScalarField H1(heatTransferModelIter().first()->K());
        const volScalarField H2(heatTransferModelIter().second()->K());
        const dimensionedScalar HSmall("small", heatTransferModel::dimK, small);

        volScalarField& Tf = *this->Tf_[pair];

        for (label i = 0; i < nInterfaceCorrectors_; ++ i)
        {
            volScalarField mDotL
            (
                IOobject
                (
                    "mDotL",
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
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
                dimensionedScalar(mDotL.dimensions()/dimTemperature, 0)
            );

            // Add latent heats from forward and backward models
            if (this->interfaceCompositionModels_.found(key12))
            {
                this->interfaceCompositionModels_[key12]->addMDotL
                (
                    massTransferModels_[pair].first()->K(),
                    Tf,
                    mDotL,
                    mDotLPrime
                );
            }
            if (this->interfaceCompositionModels_.found(key21))
            {
                this->interfaceCompositionModels_[key21]->addMDotL
                (
                    massTransferModels_[pair].second()->K(),
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
