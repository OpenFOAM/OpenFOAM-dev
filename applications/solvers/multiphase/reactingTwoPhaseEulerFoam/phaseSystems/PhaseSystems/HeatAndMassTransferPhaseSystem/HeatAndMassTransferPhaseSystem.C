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

#include "HeatAndMassTransferPhaseSystem.H"

#include "BlendedInterfacialModel.H"
#include "heatTransferModel.H"
#include "massTransferModel.H"
#include "interfaceCompositionModel.H"

#include "HashPtrTable.H"

#include "fvcDiv.H"
#include "fvmSup.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::
HeatAndMassTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generatePairsAndSubModels
    (
        "heatTransfer",
        heatTransferModels_
    );

    this->generatePairsAndSubModels
    (
        "massTransfer",
        massTransferModels_
    );

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

        // Initialy assume no mass transfer

        dmdt_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("dmdt", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimDensity/dimTime, 0)
            )
        );

        dmdtExplicit_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("dmdtExplicit", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar("zero", dimDensity/dimTime, 0)
            )
        );

        volScalarField H1(heatTransferModels_[pair][pair.first()]->K());
        volScalarField H2(heatTransferModels_[pair][pair.second()]->K());

        Tf_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Tf", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                (
                    H1*pair.phase1().thermo().T()
                  + H2*pair.phase2().thermo().T()
                )
               /max
                (
                    H1 + H2,
                    dimensionedScalar("small", heatTransferModel::dimK, SMALL)
                ),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        Tf_[pair]->correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::
~HeatAndMassTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    const scalar dmdtSign(Pair<word>::compare(dmdt_.find(key).key(), key));

    return dmdtSign**dmdt_[key];
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::momentumTransfer() const
{
    autoPtr<phaseSystem::momentumTransferTable>
        eqnsPtr(BasePhaseSystem::momentumTransfer());

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    // Source term due to mass trasfer
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

        const volScalarField dmdt(this->dmdt(pair));
        const volScalarField dmdt12(posPart(dmdt));
        const volScalarField dmdt21(negPart(dmdt));

        *eqns[pair.phase1().name()] += fvm::Sp(dmdt21, U1) - dmdt21*U2;
        *eqns[pair.phase2().name()] += dmdt12*U1 - fvm::Sp(dmdt12, U2);
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr
    (
        new phaseSystem::heatTransferTable()
    );

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    forAllConstIter
    (
        phaseSystem::phaseModelTable,
        this->phaseModels_,
        phaseModelIter
    )
    {
        const phaseModel& phase(phaseModelIter());

        eqns.insert
        (
            phase.name(),
            new fvScalarMatrix(phase.thermo().he(), dimEnergy/dimTime)
        );
    }

    // Heat transfer with the interface
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phasePair& pair
        (
            this->phasePairs_[heatTransferModelIter.key()]
        );

        const phaseModel* phase = &pair.phase1();
        const phaseModel* otherPhase = &pair.phase2();

        const volScalarField& Tf(*Tf_[pair]);

        const volScalarField K1
        (
            heatTransferModelIter()[pair.first()]->K()
        );
        const volScalarField K2
        (
            heatTransferModelIter()[pair.second()]->K()
        );
        const volScalarField KEff
        (
            K1*K2
           /max
            (
                K1 + K2,
                dimensionedScalar("small", heatTransferModel::dimK, SMALL)
            )
        );

        const volScalarField* K = &K1;
        const volScalarField* otherK = &K2;

        forAllConstIter(phasePair, pair, iter)
        {
            const volScalarField& he(phase->thermo().he());
            volScalarField Cpv(phase->thermo().Cpv());

            *eqns[phase->name()] +=
                (*K)*(Tf - phase->thermo().T())
              + KEff/Cpv*he - fvm::Sp(KEff/Cpv, he);

            Swap(phase, otherPhase);
            Swap(K, otherK);
        }
    }

    // Source term due to mass trasfer
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

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        const volScalarField& he1(phase1.thermo().he());
        const volScalarField& he2(phase2.thermo().he());

        const volScalarField& K1(phase1.K());
        const volScalarField& K2(phase2.K());

        const volScalarField dmdt(this->dmdt(pair));
        const volScalarField dmdt12(posPart(dmdt));
        const volScalarField dmdt21(negPart(dmdt));
        const volScalarField& Tf(*Tf_[pair]);

        *eqns[phase1.name()] +=
            fvm::Sp(dmdt21, he1) + dmdt21*K1
          - dmdt21*(phase2.thermo().he(phase2.thermo().p(), Tf) + K2);

        *eqns[phase2.name()] +=
            dmdt12*(phase1.thermo().he(phase1.thermo().p(), Tf) + K1)
          - fvm::Sp(dmdt12, he2) - dmdt12*K2;
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::massTransferTable>
Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::massTransfer() const
{
    // Create a mass transfer matrix for each species of each phase
    autoPtr<phaseSystem::massTransferTable> eqnsPtr
    (
        new phaseSystem::massTransferTable()
    );

    phaseSystem::massTransferTable& eqns = eqnsPtr();

    forAllConstIter
    (
        phaseSystem::phaseModelTable,
        this->phaseModels_,
        phaseModelIter
    )
    {
        const phaseModel& phase(phaseModelIter());

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

        *dmdt_[pair] =
            *dmdtExplicit_[pair];

        *dmdtExplicit_[pair] =
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

        const volScalarField& Tf(*Tf_[key]);

        volScalarField& dmdtExplicit(*dmdtExplicit_[key]);
        volScalarField& dmdt(*dmdt_[key]);

        scalar dmdtSign(Pair<word>::compare(dmdt_.find(key).key(), key));

        const volScalarField K
        (
            massTransferModels_[key][phase.name()]->K()
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
void Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::correctThermo()
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

        volScalarField H1(heatTransferModels_[pair][pair.first()]->K());
        volScalarField H2(heatTransferModels_[pair][pair.second()]->K());
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

        volScalarField& Tf = *Tf_[pair];

        // Add latent heats from forward and backward models
        if (interfaceCompositionModels_.found(key12))
        {
            interfaceCompositionModels_[key12]->addMDotL
            (
                massTransferModels_[pair][pair.first()]->K(),
                Tf,
                mDotL,
                mDotLPrime
            );
        }
        if (interfaceCompositionModels_.found(key21))
        {
            interfaceCompositionModels_[key21]->addMDotL
            (
                massTransferModels_[pair][pair.second()]->K(),
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

        // Update the interface compositions
        if (interfaceCompositionModels_.found(key12))
        {
            interfaceCompositionModels_[key12]->update(Tf);
        }
        if (interfaceCompositionModels_.found(key21))
        {
            interfaceCompositionModels_[key21]->update(Tf);
        }

        Tf.correctBoundaryConditions();

        Info<< "Tf." << pair.name()
            << ": min = " << min(Tf.internalField())
            << ", mean = " << average(Tf.internalField())
            << ", max = " << max(Tf.internalField())
            << endl;
    }
}


template<class BasePhaseSystem>
bool Foam::HeatAndMassTransferPhaseSystem<BasePhaseSystem>::read()
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
