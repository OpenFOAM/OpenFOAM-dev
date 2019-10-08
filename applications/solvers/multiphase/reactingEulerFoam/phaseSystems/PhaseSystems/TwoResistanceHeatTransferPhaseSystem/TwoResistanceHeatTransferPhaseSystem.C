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

#include "TwoResistanceHeatTransferPhaseSystem.H"
#include "BlendedInterfacialModel.H"
#include "heatTransferModel.H"
#include "fvmSup.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::addDmdtHefs
(
    const phaseSystem::dmdtfTable& dmdtfs,
    phaseSystem::heatTransferTable& eqns
) const
{
    // In the presence of thermally-coupled mass transfer, the relation between
    // heat transfers across the interface between phases 1 and 2 is:
    //
    //                         Q1 + Q2 = mDot*L
    //     H1*(Tf - T1) + H2*(Tf - T1) = mDot*L
    //
    // Where Q1 and Q2 are the net transfer into phases 1 and 2 respectively,
    // H1 and H2 are the heat transfer coefficients on either side, Tf is the
    // temperature at the interface, mDot is the mass transfer rate from phase
    // 2 to phase 1, and L is the latent heat of phase 2 minus phase 1.
    //
    // Rearranging for Tf:
    //
    //    Tf = (H1*T1 + H2*T2 + mDot*L)/(H1 + H2)
    //
    // And for Q1 and Q2:
    //
    //    Q1 = HEff*(T2 - T1) + H1/(H1 + H2)*mDot*L
    //    Q2 = HEff*(T1 - T2) + H2/(H1 + H2)*mDot*L
    //
    // The HEff terms are those implemented in the heatTransfer method. The
    // latent heat mDot*L terms are added below (as well as other standard
    // terms representing sensible energy and kinetic energy differences).

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
        // MomentumTransferPhaseSystem::addDmdtU for an explanation of the
        // fvm::Sp terms below.

        if (heatTransferModels_.found(key))
        {
            // Assume a thermally-coupled mass transfer process. Calculate
            // heat transfers by evaluating properties at the interface
            // temperature...

            // Transfer of energy from the interface into the bulk
            const volScalarField& Tf(*Tf_[key]);
            const volScalarField hef1(thermo1.he(thermo1.p(), Tf));
            const volScalarField hef2(thermo2.he(thermo2.p(), Tf));
            *eqns[phase1.name()] += dmdtf*hef1 - fvm::Sp(dmdtf, he1);
            *eqns[phase2.name()] -= dmdtf*hef2 - fvm::Sp(dmdtf, he2);

            // Latent heat contribution
            const volScalarField L(hef2 + thermo2.hc() - hef1 - thermo1.hc());
            const volScalarField H1(heatTransferModels_[key].first()->K());
            const volScalarField H2(heatTransferModels_[key].second()->K());
            const volScalarField H1Fac(H1/(H1 + H2));
            *eqns[phase1.name()] += H1Fac*dmdtf*L;
            *eqns[phase2.name()] += (1 - H1Fac)*dmdtf*L;

            // Transfer of kinetic energy
            *eqns[phase1.name()] += dmdtf21*(K2 - K1);
            *eqns[phase2.name()] -= dmdtf12*(K1 - K2);
        }
        else
        {
            // In the absence of a heat transfer model, assume a
            // non-thermally-coupled mass transfer process; i.e., no phase
            // change and therefore no interface state or latent heat...

            // Transfer of energy from bulk to bulk
            *eqns[phase1.name()] += dmdtf21*he2 - fvm::Sp(dmdtf21, he1);
            *eqns[phase2.name()] -= dmdtf12*he1 - fvm::Sp(dmdtf12, he2);

            // Transfer of kinetic energy
            *eqns[phase1.name()] += dmdtf21*(K2 - K1);
            *eqns[phase2.name()] -= dmdtf12*(K1 - K2);
        }
    }
}


template<class BasePhaseSystem>
void Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::addDmidtHef
(
    const phaseSystem::dmidtfTable& dmidtfs,
    phaseSystem::heatTransferTable& eqns
) const
{
    // See the addDmdtHe method for a description of the latent heat terms

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
        // MomentumTransferPhaseSystem::addDmdtU for an explanation of the
        // fvm::Sp terms below.

        if (heatTransferModels_.found(key))
        {
            // Assume a thermally-coupled mass transfer process. Calculate
            // heat transfers by evaluating properties at the interface
            // temperature...

            // Transfer coefficients
            const volScalarField H1(heatTransferModels_[key].first()->K());
            const volScalarField H2(heatTransferModels_[key].second()->K());
            const volScalarField H1Fac(H1/(H1 + H2));

            // Interface properties
            const volScalarField& Tf(*Tf_[key]);
            const volScalarField hef1(thermo1.he(thermo1.p(), Tf));
            const volScalarField hef2(thermo2.he(thermo2.p(), Tf));
            const volScalarField hc1(thermo1.hc());
            const volScalarField hc2(thermo2.hc());

            // Loop the species
            forAllConstIter
            (
                HashPtrTable<volScalarField>,
                *dmidtfIter(),
                dmidtfJter
            )
            {
                const word& member = dmidtfJter.key();

                const volScalarField dmidtf
                (
                    Pair<word>::compare(pair, key)**dmidtfJter()
                );
                const volScalarField dmidtf21(posPart(dmidtf));
                const volScalarField dmidtf12(negPart(dmidtf));

                // Create the interface energies for the transferring specie
                volScalarField hefi1(hef1), hci1(hc1);
                if (isA<rhoReactionThermo>(thermo1))
                {
                    const basicSpecieMixture& composition1 =
                        refCast<const rhoReactionThermo>(thermo1).composition();
                    hefi1 =
                        composition1.HE
                        (
                            composition1.species()[member],
                            thermo1.p(),
                            Tf
                        );
                    hci1 =
                        dimensionedScalar
                        (
                            dimEnergy/dimMass,
                            composition1.Hc(composition1.species()[member])
                        );
                }
                volScalarField hefi2(hef2), hci2(hc2);
                if (isA<rhoReactionThermo>(thermo2))
                {
                    const basicSpecieMixture& composition2 =
                        refCast<const rhoReactionThermo>(thermo2).composition();
                    hefi2 =
                        composition2.HE
                        (
                            composition2.species()[member],
                            thermo2.p(),
                            Tf
                        );
                    hci2 =
                        dimensionedScalar
                        (
                            dimEnergy/dimMass,
                            composition2.Hc(composition2.species()[member])
                        );
                }

                // Create the latent heat for the transferring specie
                const volScalarField Li(hefi2 + hci2 - hefi1 - hci1);

                // Transfer of energy from the interface into the bulk
                *eqns[phase1.name()] += dmidtf*hefi1 - fvm::Sp(dmidtf, he1);
                *eqns[phase2.name()] -= dmidtf*hefi2 - fvm::Sp(dmidtf, he2);

                // Latent heat contribution
                *eqns[phase1.name()] += H1Fac*dmidtf*Li;
                *eqns[phase2.name()] += (1 - H1Fac)*dmidtf*Li;

                // Transfer of kinetic energy
                *eqns[phase1.name()] += dmidtf21*(K2 - K1);
                *eqns[phase2.name()] -= dmidtf12*(K1 - K2);
            }
        }
        else
        {
            // In the absence of a heat transfer model, assume a
            // non-thermally-coupled mass transfer process; i.e., no phase
            // change and therefore no interface state or latent heat...

            // Loop the species
            forAllConstIter
            (
                HashPtrTable<volScalarField>,
                *dmidtfIter(),
                dmidtfJter
            )
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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
TwoResistanceHeatTransferPhaseSystem
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

    // Check that models have been specified on both sides of the interfaces
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[heatTransferModelIter.key()];

        if (!heatTransferModels_[pair].first().valid())
        {
            FatalErrorInFunction
                << "A heat transfer model for the " << pair.phase1().name()
                << " side of the " << pair << " pair is not specified"
                << exit(FatalError);
        }
        if (!heatTransferModels_[pair].second().valid())
        {
            FatalErrorInFunction
                << "A heat transfer model for the " << pair.phase2().name()
                << " side of the " << pair << " pair is not specified"
                << exit(FatalError);
        }
    }

    // Calculate initial Tf-s as the mean between the two phases
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[heatTransferModelIter.key()];

        const volScalarField& T1(pair.phase1().thermo().T());
        const volScalarField& T2(pair.phase2().thermo().T());

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
                (T1 + T2)/2
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
~TwoResistanceHeatTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
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

    // In the absence of mass transfer, the relation between heat transfers
    // across the interface between phases 1 and 2 is:
    //
    //                         Q1 + Q2 = 0
    //     H1*(Tf - T1) + H2*(Tf - T1) = 0
    //
    // Where Q1 and Q2 are the net transfer into phases 1 and 2 respectively,
    // H1 and H2 are the heat transfer coefficients on either side, and Tf is
    // the temperature at the interface.
    //
    // Rearranging for Tf:
    //
    //     Tf = (H1*T1 + H2*T2)/(H1 + H2)
    //
    // Which gives the following for Q1 and Q2:
    //
    //     Q1 = H1*H2/(H1 + H2)*(T2 - T1)
    //     Q2 = H1*H2/(H1 + H2)*(T1 - T2)
    //
    // The transfer is therefore the same as for the single-resistance system
    // with a single effective heat transfer coefficient:
    //
    //     HEff = H1*H2/(H1 + H2)
    //
    // This is implemented below

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

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        const volScalarField& he1(phase1.thermo().he());
        const volScalarField& he2(phase2.thermo().he());

        const volScalarField Cpv1(phase1.thermo().Cpv());
        const volScalarField Cpv2(phase2.thermo().Cpv());

        // Heat transfer between phases in the absence of mass transfer
        const volScalarField H1(heatTransferModelIter().first()->K());
        const volScalarField H2(heatTransferModelIter().second()->K());
        const volScalarField HEff(H1*H2/(H1 + H2));
        *eqns[phase1.name()] +=
            HEff*(phase2.thermo().T() - phase1.thermo().T())
          + H1/Cpv1*he1 - fvm::Sp(H1/Cpv1, he1);
        *eqns[phase2.name()] +=
            HEff*(phase1.thermo().T() - phase2.thermo().T())
          + H2/Cpv2*he2 - fvm::Sp(H2/Cpv2, he2);
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
correctEnergyTransport()
{
    BasePhaseSystem::correctEnergyTransport();

    correctInterfaceThermo();
}


template<class BasePhaseSystem>
void Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
correctInterfaceThermo()
{
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[heatTransferModelIter.key()];

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        const volScalarField& T1(phase1.thermo().T());
        const volScalarField& T2(phase2.thermo().T());

        const volScalarField H1(heatTransferModels_[pair].first()->K());
        const volScalarField H2(heatTransferModels_[pair].second()->K());
        const dimensionedScalar HSmall(H1.dimensions(), small);

        volScalarField& Tf = *Tf_[pair];

        Tf = (H1*T1 + H2*T2)/max(H1 + H2, HSmall);

        Info<< "Tf." << pair.name()
            << ": min = " << min(Tf.primitiveField())
            << ", mean = " << average(Tf.primitiveField())
            << ", max = " << max(Tf.primitiveField())
            << endl;
    }
}


template<class BasePhaseSystem>
bool Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::read()
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
