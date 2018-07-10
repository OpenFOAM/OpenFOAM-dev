/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "addToRunTimeSelectionTable.H"

#include "phaseSystem.H"
#include "multiphaseSystem.H"
#include "MomentumTransferPhaseSystem.H"
#include "OneResistanceHeatTransferPhaseSystem.H"
#include "TwoResistanceHeatTransferPhaseSystem.H"
#include "PhaseTransferPhaseSystem.H"
#include "PopulationBalancePhaseSystem.H"
#include "InterfaceCompositionPhaseChangePhaseSystem.H"
#include "ThermalPhaseChangePhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        PhaseTransferPhaseSystem
        <
            OneResistanceHeatTransferPhaseSystem
            <
                MomentumTransferPhaseSystem<multiphaseSystem>
            >
        >
        basicMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        multiphaseSystem,
        basicMultiphaseSystem,
        dictionary,
        basicMultiphaseSystem
    );

    typedef
        InterfaceCompositionPhaseChangePhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                TwoResistanceHeatTransferPhaseSystem
                <
                    MomentumTransferPhaseSystem<multiphaseSystem>
                >
            >
        >
        interfaceCompositionPhaseChangeMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        multiphaseSystem,
        interfaceCompositionPhaseChangeMultiphaseSystem,
        dictionary,
        interfaceCompositionPhaseChangeMultiphaseSystem
    );

    typedef
        ThermalPhaseChangePhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                TwoResistanceHeatTransferPhaseSystem
                <
                    MomentumTransferPhaseSystem<multiphaseSystem>
                >
            >
        >
        thermalPhaseChangeMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        multiphaseSystem,
        thermalPhaseChangeMultiphaseSystem,
        dictionary,
        thermalPhaseChangeMultiphaseSystem
    );

    typedef
        PopulationBalancePhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                OneResistanceHeatTransferPhaseSystem
                <
                    MomentumTransferPhaseSystem<multiphaseSystem>
                >
            >
        >
        populationBalanceMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        multiphaseSystem,
        populationBalanceMultiphaseSystem,
        dictionary,
        populationBalanceMultiphaseSystem
    );

    typedef
        ThermalPhaseChangePhaseSystem
        <
            PopulationBalancePhaseSystem
            <
                PhaseTransferPhaseSystem
                <
                    TwoResistanceHeatTransferPhaseSystem
                    <
                        MomentumTransferPhaseSystem<multiphaseSystem>
                    >
                >
            >
        >
        thermalPhaseChangePopulationBalanceMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        multiphaseSystem,
        thermalPhaseChangePopulationBalanceMultiphaseSystem,
        dictionary,
        thermalPhaseChangePopulationBalanceMultiphaseSystem
    );
}


// ************************************************************************* //
