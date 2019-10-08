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

#include "OneResistanceHeatTransferPhaseSystem.H"
#include "BlendedInterfacialModel.H"
#include "heatTransferModel.H"
#include "fvmSup.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::addDmdtHefs
(
    const phaseSystem::dmdtfTable& dmdtfs,
    phaseSystem::heatTransferTable& eqns
) const
{
    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs, dmdtfIter)
    {
        const phasePairKey& key = dmdtfIter.key();
        const phasePair& pair(this->phasePairs_[key]);
        const volScalarField dmdtf(Pair<word>::compare(pair, key)**dmdtfIter());
        const volScalarField dmdtf21(posPart(dmdtf));
        const volScalarField dmdtf12(negPart(dmdtf));

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();
        const rhoThermo& thermo1 = phase1.thermo();
        const rhoThermo& thermo2 = phase2.thermo();
        const volScalarField& he1(thermo1.he());
        const volScalarField& he2(thermo2.he());
        const volScalarField K1(phase1.K());
        const volScalarField K2(phase2.K());

        // Note that the phase EEqn contains a continuity error term. See
        // MomentumTransferPhaseSystem::addDmdtUfs for an explanation of the
        // fvm::Sp terms below.

        // Transfer of energy from bulk to bulk
        *eqns[phase1.name()] += dmdtf21*he2 - fvm::Sp(dmdtf21, he1);
        *eqns[phase2.name()] -= dmdtf12*he1 - fvm::Sp(dmdtf12, he2);

        // Transfer of kinetic energy
        *eqns[phase1.name()] += dmdtf21*(K2 - K1);
        *eqns[phase2.name()] -= dmdtf12*(K1 - K2);
    }
}


template<class BasePhaseSystem>
void Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::addDmidtHef
(
    const phaseSystem::dmidtfTable& dmidtfs,
    phaseSystem::heatTransferTable& eqns
) const
{
    forAllConstIter(phaseSystem::dmidtfTable, dmidtfs, dmidtfIter)
    {
        const phasePairKey& key = dmidtfIter.key();
        const phasePair& pair(this->phasePairs_[key]);

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();
        const rhoThermo& thermo1 = phase1.thermo();
        const rhoThermo& thermo2 = phase2.thermo();
        const volScalarField& he1(thermo1.he());
        const volScalarField& he2(thermo2.he());
        const volScalarField K1(phase1.K());
        const volScalarField K2(phase2.K());

        // Note that the phase EEqn contains a continuity error term. See
        // MomentumTransferPhaseSystem::addDmdtUfs for an explanation of the
        // fvm::Sp terms below.

        forAllConstIter(HashPtrTable<volScalarField>, *dmidtfIter(), dmidtfJter)
        {
            const word& member = dmidtfJter.key();

            const volScalarField dmidtf
            (
                Pair<word>::compare(pair, key)**dmidtfJter()
            );
            const volScalarField dmidtf21(posPart(dmidtf));
            const volScalarField dmidtf12(negPart(dmidtf));

            // Create the energies for the transferring specie
            volScalarField hei1(he1);
            if (isA<rhoReactionThermo>(thermo1))
            {
                const basicSpecieMixture& composition1 =
                    refCast<const rhoReactionThermo>(thermo1).composition();
                hei1 =
                    composition1.HE
                    (
                        composition1.species()[member],
                        thermo1.p(),
                        thermo1.T()
                    );
            }
            volScalarField hei2(he2);
            if (isA<rhoReactionThermo>(thermo2))
            {
                const basicSpecieMixture& composition2 =
                    refCast<const rhoReactionThermo>(thermo2).composition();
                hei2 =
                    composition2.HE
                    (
                        composition2.species()[member],
                        thermo2.p(),
                        thermo2.T()
                    );
            }

            // Transfer of energy from bulk to bulk
            *eqns[phase1.name()] += dmidtf21*hei2 - fvm::Sp(dmidtf21, he1);
            *eqns[phase2.name()] -= dmidtf12*hei1 - fvm::Sp(dmidtf12, he2);

            // Transfer of kinetic energy
            *eqns[phase1.name()] += dmidtf21*(K2 - K1);
            *eqns[phase2.name()] -= dmidtf12*(K1 - K2);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
OneResistanceHeatTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generatePairsAndSubModels
    (
        "heatTransfer",
        heatTransferModels_,
        false
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
~OneResistanceHeatTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr
    (
        new phaseSystem::heatTransferTable()
    );

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        eqns.insert
        (
            phase.name(),
            new fvScalarMatrix(phase.thermo().he(), dimEnergy/dimTime)
        );
    }

    // Heat transfer across the interface
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const volScalarField K(heatTransferModelIter()->K());

        const phasePair& pair(this->phasePairs_[heatTransferModelIter.key()]);

        forAllConstIter(phasePair, pair, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const volScalarField& he(phase.thermo().he());
            volScalarField Cpv(phase.thermo().Cpv());

            *eqns[phase.name()] +=
                K*(otherPhase.thermo().T() - phase.thermo().T() + he/Cpv)
              - fvm::Sp(K/Cpv, he);
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
bool Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::read()
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
