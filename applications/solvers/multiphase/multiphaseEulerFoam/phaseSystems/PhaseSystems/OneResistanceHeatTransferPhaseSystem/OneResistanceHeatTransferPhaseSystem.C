/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::OneResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
OneResistanceHeatTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    HeatTransferPhaseSystem<BasePhaseSystem>(mesh)
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

            const volScalarField& he = phase.thermo().he();
            const volScalarField Cpv(phase.thermo().Cpv());

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
    if (HeatTransferPhaseSystem<BasePhaseSystem>::read())
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
