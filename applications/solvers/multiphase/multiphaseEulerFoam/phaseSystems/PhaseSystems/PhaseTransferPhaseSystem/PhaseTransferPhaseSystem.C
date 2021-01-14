/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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

#include "PhaseTransferPhaseSystem.H"
#include "phaseTransferModel.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::dmdtfTable>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::totalDmdtfs() const
{
    autoPtr<phaseSystem::dmdtfTable> totalDmdtfsPtr
    (
        new phaseSystem::dmdtfTable
    );
    phaseSystem::dmdtfTable& totalDmdtfs = totalDmdtfsPtr();

    forAllConstIter
    (
        phaseTransferModelTable,
        phaseTransferModels_,
        phaseTransferModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[phaseTransferModelIter.key()];

        totalDmdtfs.insert(pair, phaseSystem::dmdtf(pair).ptr());

        if (phaseTransferModelIter()->mixture())
        {
            *totalDmdtfs[pair] += *dmdtfs_[pair];
        }

        forAllConstIter
        (
            HashPtrTable<volScalarField>,
            *dmidtfs_[pair],
            dmidtfIter
        )
        {
            *totalDmdtfs[pair] += *dmidtfIter();
        }
    }

    return totalDmdtfsPtr;
}


// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::addDmdtYfs
(
    const phaseSystem::dmdtfTable& dmdtfs,
    phaseSystem::specieTransferTable& eqns
) const
{
    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs, dmdtfIter)
    {
        const phasePairKey& key = dmdtfIter.key();
        const phasePair& pair(this->phasePairs_[key]);

        const volScalarField dmdtf(Pair<word>::compare(pair, key)**dmdtfIter());
        const volScalarField dmdtf12(negPart(dmdtf));
        const volScalarField dmdtf21(posPart(dmdtf));

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        forAll(phase1.Y(), Yi1)
        {
            const volScalarField& Y1 = phase1.Y()[Yi1];
            const volScalarField& Y2 = phase2.Y(Y1.member());

            *eqns[Y1.name()] += dmdtf21*Y2 + fvm::Sp(dmdtf12, Y1);
            *eqns[Y2.name()] -= dmdtf12*Y1 + fvm::Sp(dmdtf21, Y2);
        }
    }
}


template<class BasePhaseSystem>
void Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::addDmidtYf
(
    const phaseSystem::dmidtfTable& dmidtfs,
    phaseSystem::specieTransferTable& eqns
) const
{
    forAllConstIter(phaseSystem::dmidtfTable, dmidtfs, dmidtfIter)
    {
        const phasePairKey& key = dmidtfIter.key();
        const phasePair& pair(this->phasePairs_[key]);

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        forAllConstIter(HashPtrTable<volScalarField>, *dmidtfIter(), dmidtfJter)
        {
            const word& member = dmidtfJter.key();

            const volScalarField dmidtf
            (
                Pair<word>::compare(pair, key)**dmidtfJter()
            );

            if (!phase1.pure())
            {
                const volScalarField& Y1 = phase1.Y(member);
                *eqns[Y1.name()] += dmidtf;
            }

            if (!phase2.pure())
            {
                const volScalarField& Y2 = phase2.Y(member);
                *eqns[Y2.name()] -= dmidtf;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::PhaseTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generatePairsAndSubModels
    (
        "phaseTransfer",
        phaseTransferModels_,
        false
    );

    forAllConstIter
    (
        phaseTransferModelTable,
        phaseTransferModels_,
        phaseTransferModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[phaseTransferModelIter.key()];

        this->template validateMassTransfer<phaseTransferModel>(pair);

        if (phaseTransferModelIter()->mixture())
        {
            dmdtfs_.insert
            (
                pair,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "phaseTransfer:dmdtf",
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

        dmidtfs_.insert(pair, new HashPtrTable<volScalarField>());

        const hashedWordList species(phaseTransferModelIter()->species());

        forAllConstIter(hashedWordList, species, specieIter)
        {
            const word& specie = *specieIter;

            dmidtfs_[pair]->insert
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
                                "phaseTransfer:dmidtf",
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


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::~PhaseTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::dmdtf
(
    const phasePairKey& key
) const
{
    tmp<volScalarField> tDmdtf = BasePhaseSystem::dmdtf(key);

    if (phaseTransferModels_.found(key))
    {
        const label dmdtfSign(Pair<word>::compare(this->phasePairs_[key], key));

        if (phaseTransferModels_[key]->mixture())
        {
            tDmdtf.ref() += dmdtfSign**dmdtfs_[key];
        }

        forAllConstIter
        (
            HashPtrTable<volScalarField>,
            *dmidtfs_[key],
            dmidtfIter
        )
        {
            tDmdtf.ref() += dmdtfSign**dmidtfIter();
        }
    }

    return tDmdtf;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    autoPtr<phaseSystem::dmdtfTable> totalDmdtfsPtr = this->totalDmdtfs();
    const phaseSystem::dmdtfTable& totalDmdtfs = totalDmdtfsPtr();

    forAllConstIter(phaseSystem::dmdtfTable, totalDmdtfs, totalDmdtfIter)
    {
        const phasePair& pair = this->phasePairs_[totalDmdtfIter.key()];

        addField(pair.phase1(), "dmdt", *totalDmdtfIter(), dmdts);
        addField(pair.phase2(), "dmdt", - *totalDmdtfIter(), dmdts);
    }

    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::momentumTransfer()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransfer();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(totalDmdtfs(), eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::momentumTransferf()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransferf();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(totalDmdtfs(), eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    this->addDmdtHefs(dmdtfs_, eqns);
    this->addDmidtHefs(dmidtfs_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::specieTransferTable>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::specieTransfer() const
{
    autoPtr<phaseSystem::specieTransferTable> eqnsPtr
    (
        new phaseSystem::specieTransferTable()
    );

    phaseSystem::specieTransferTable& eqns = eqnsPtr();

    // Create a mass transfer matrix for each species of each phase
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

    this->addDmdtYfs(dmdtfs_, eqns);
    this->addDmidtYf(dmidtfs_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::correct()
{
    BasePhaseSystem::correct();

    // Reset all the mass transfer rates to zero
    forAllConstIter
    (
        phaseTransferModelTable,
        phaseTransferModels_,
        phaseTransferModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[phaseTransferModelIter.key()];

        if (phaseTransferModelIter()->mixture())
        {
            *dmdtfs_[pair] = Zero;
        }

        const hashedWordList species(phaseTransferModelIter()->species());

        forAllConstIter(hashedWordList, species, specieIter)
        {
            const word& specie = *specieIter;

            *(*dmidtfs_[pair])[specie] = Zero;
        }
    }

    // Evaluate the models and sum the results into the mass transfer tables
    forAllIter
    (
        phaseTransferModelTable,
        phaseTransferModels_,
        phaseTransferModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[phaseTransferModelIter.key()];

        if (phaseTransferModelIter()->mixture())
        {
            *dmdtfs_[pair] += phaseTransferModelIter()->dmdtf();
        }

        const HashPtrTable<volScalarField> dmidtf
        (
            phaseTransferModelIter()->dmidtf()
        );

        forAllConstIter(HashPtrTable<volScalarField>, dmidtf, dmidtfIter)
        {
            *(*dmidtfs_[pair])[dmidtfIter.key()] += *dmidtfIter();
        }
    }
}


template<class BasePhaseSystem>
bool Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::read()
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
