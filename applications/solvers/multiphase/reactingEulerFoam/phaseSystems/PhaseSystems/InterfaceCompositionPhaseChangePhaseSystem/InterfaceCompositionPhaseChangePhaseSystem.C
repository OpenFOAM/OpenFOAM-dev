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
Foam::autoPtr<Foam::phaseSystem::dmdtfTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
dmdtfs() const
{
    autoPtr<phaseSystem::dmdtfTable> dmdtPtr(new phaseSystem::dmdtfTable);
    phaseSystem::dmdtfTable& dmdtfs = dmdtPtr();

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
                specieIter
            )
            {
                const word& specie = *specieIter;

                tmp<volScalarField> dmidtf
                (
                    (pairIter.index() == 0 ? +1 : -1)
                   *(
                        *(*dmidtfSus_[pair])[specie]
                      + *(*dmidtfSps_[pair])[specie]*phase.Y(specie)
                    )
                );

                if (dmdtfs.found(pair))
                {
                    *dmdtfs[pair] += dmidtf;
                }
                else
                {
                    dmdtfs.insert(pair, dmidtf.ptr());
                }
            }
        }
    }

    return dmdtPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::dmidtfTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
dmidtfs() const
{
    autoPtr<phaseSystem::dmidtfTable> dmidtfsPtr(new phaseSystem::dmidtfTable);
    phaseSystem::dmidtfTable& dmidtfs = dmidtfsPtr();

    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];

        if (!dmidtfs.found(pair))
        {
            dmidtfs.insert(pair, new HashPtrTable<volScalarField>());
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
                specieIter
            )
            {
                const word& specie = *specieIter;


                tmp<volScalarField> dmidtf
                (
                    (pairIter.index() == 0 ? +1 : -1)
                   *(
                        *(*dmidtfSus_[pair])[specie]
                      + *(*dmidtfSps_[pair])[specie]*phase.Y(specie)
                    )
                );

                if (dmidtfs[pair]->found(specie))
                {
                    *(*dmidtfs[pair])[specie] += dmidtf;
                }
                else
                {
                    dmidtfs[pair]->insert(specie, dmidtf.ptr());
                }
            }
        }
    }

    return dmidtfsPtr;
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

        dmidtfSus_.insert(pair, new HashPtrTable<volScalarField>());
        dmidtfSps_.insert(pair, new HashPtrTable<volScalarField>());

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
                specieIter
            )
            {
                const word& specie = *specieIter;

                dmidtfSus_[pair]->insert
                (
                    specie,
                    new volScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                IOobject::groupName
                                (
                                    "interfaceCompositionPhaseChange:dmidtfSu",
                                    specie
                                ),
                                pair.name()
                            ),
                            this->mesh().time().timeName(),
                            this->mesh()
                        ),
                        this->mesh(),
                        dimensionedScalar(dimDensity/dimTime, 0)
                    )
                );

                dmidtfSps_[pair]->insert
                (
                    specie,
                    new volScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                IOobject::groupName
                                (
                                    "interfaceCompositionPhaseChange:dmidtfSp",
                                    specie
                                ),
                                pair.name()
                            ),
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
~InterfaceCompositionPhaseChangePhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::dmdtf
(
    const phasePairKey& key
) const
{
    tmp<volScalarField> tDmdtf = BasePhaseSystem::dmdtf(key);

    if (interfaceCompositionModels_.found(key))
    {
        const phasePair& pair = this->phasePairs_[key];

        forAllConstIter(phasePair, pair, pairIter)
        {
            const autoPtr<interfaceCompositionModel>& compositionModelPtr =
                interfaceCompositionModels_[key][pairIter.index()];

            if (!compositionModelPtr.valid()) continue;

            const interfaceCompositionModel& compositionModel =
                compositionModelPtr();

            const phaseModel& phase = *pairIter;

            forAllConstIter
            (
                hashedWordList,
                compositionModel.species(),
                specieIter
            )
            {
                const word& specie = *specieIter;

                tDmdtf.ref() +=
                    Pair<word>::compare(pair, key)
                   *(
                        *(*dmidtfSus_[pair])[specie]
                      - *(*dmidtfSps_[pair])[specie]*phase.Y(specie)
                    );
            }
        }
    }

    return tDmdtf;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    autoPtr<phaseSystem::dmdtfTable> dmdtfsPtr = this->dmdtfs();
    const phaseSystem::dmdtfTable& dmdtfs = dmdtfsPtr();

    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs, dmdtfsIter)
    {
        const phasePair& pair = this->phasePairs_[dmdtfsIter.key()];
        const phaseModel& phase = pair.phase1();
        const phaseModel& otherPhase = pair.phase2();

        addField(phase, "dmdt", *dmdtfsIter(), dmdts);
        addField(otherPhase, "dmdt", - *dmdtfsIter(), dmdts);
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

    this->addDmdtUfs(dmdtfs(), eqns);

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

    this->addDmdtUfs(dmdtfs(), eqns);

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

    this->addDmidtHef(dmidtfs(), eqns);

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
    this->addDmidtYf(dmidtfs(), eqns);
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
                specieIter
            )
            {
                const word& specie = *specieIter;

                // Implicit transport through this phase
                *eqns[phase.Y(specie).name()] +=
                    *(*dmidtfSus_[pair])[specie]
                  + fvm::Sp(*(*dmidtfSps_[pair])[specie], phase.Y(specie));

                // Explicit transport out of the other phase
                if (eqns.found(IOobject::groupName(specie, otherPhase.name())))
                {
                    *eqns[otherPhase.Y(specie).name()] -=
                        *(*dmidtfSus_[pair])[specie]
                      + *(*dmidtfSps_[pair])[specie]*phase.Y(specie);
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
                specieIter
            )
            {
                const word& specie = *specieIter;

                const volScalarField KD(K*compositionModel.D(specie));
                const volScalarField Yf(compositionModel.Yf(specie, Tf));

                *(*dmidtfSus_[pair])[specie] = phase.rho()*KD*Yf;
                *(*dmidtfSps_[pair])[specie] = - phase.rho()*KD;
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
            tmp<volScalarField> dmdtLf =
                volScalarField::New
                (
                    IOobject::groupName("dmdtLf", pair.name()),
                    this->mesh(),
                    dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
                );
            tmp<volScalarField> dmdtLfPrime =
                volScalarField::New
                (
                    IOobject::groupName("dmdtLfPrime", pair.name()),
                    this->mesh(),
                    dimensionedScalar(dmdtLf().dimensions()/dimTemperature, 0)
                );

            // Add latent heats from forward and backward models
            if (this->interfaceCompositionModels_[pair].first().valid())
            {
                this->interfaceCompositionModels_[pair].first()->addDmdtL
                (
                    diffusiveMassTransferModels_[pair].first()->K(),
                    Tf,
                    dmdtLf.ref(),
                    dmdtLfPrime.ref()
                );
            }
            if (this->interfaceCompositionModels_[pair].second().valid())
            {
                this->interfaceCompositionModels_[pair].second()->addDmdtL
                (
                  - diffusiveMassTransferModels_[pair].second()->K(),
                    Tf,
                    dmdtLf.ref(),
                    dmdtLfPrime.ref()
                );
            }

            // Update the interface temperature by applying one step of newton's
            // method to the interface relation
            Tf -=
                (
                    H1*(Tf - pair.phase1().thermo().T())
                  + H2*(Tf - pair.phase2().thermo().T())
                  - dmdtLf
                )
               /(
                    max(H1 + H2 - dmdtLfPrime, HSmall)
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
