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

#include "HeatTransferPhaseSystem.H"

#include "BlendedInterfacialModel.H"
#include "heatTransferModel.H"

#include "HashPtrTable.H"

#include "fvcDiv.H"
#include "fvmSup.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::HeatTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generatePairsAndSubModels("heatTransfer", heatTransferModels_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::~HeatTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
bool Foam::HeatTransferPhaseSystem<BasePhaseSystem>::transfersMass
(
    const phaseModel& phase
) const
{
    return false;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName
                (
                    "dmdt",
                    this->phasePairs_[key]->name()
                ),
                this->mesh().time().timeName(),
                this->mesh().time()
            ),
            this->mesh(),
            dimensionedScalar("zero", dimDensity/dimTime, 0)
        )
    );
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::dmdt
(
    const Foam::phaseModel& phase
) const
{
    return tmp<volScalarField>(NULL);
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::heatTransfer() const
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

    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const volScalarField K(heatTransferModelIter()->K());

        const phasePair& pair(this->phasePairs_[heatTransferModelIter.key()]);

        const phaseModel* phase = &pair.phase1();
        const phaseModel* otherPhase = &pair.phase2();

        forAllConstIter(phasePair, pair, iter)
        {
            const volScalarField& he(phase->thermo().he());
            volScalarField Cpv(phase->thermo().Cpv());

            *eqns[phase->name()] +=
                K*(otherPhase->thermo().T() - phase->thermo().T() + he/Cpv)
              - fvm::Sp(K/Cpv, he);

            Swap(phase, otherPhase);
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::massTransferTable>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::massTransfer() const
{
    autoPtr<phaseSystem::massTransferTable> eqnsPtr
    (
        new phaseSystem::massTransferTable()
    );

    return eqnsPtr;
}


template<class BasePhaseSystem>
bool Foam::HeatTransferPhaseSystem<BasePhaseSystem>::read()
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
