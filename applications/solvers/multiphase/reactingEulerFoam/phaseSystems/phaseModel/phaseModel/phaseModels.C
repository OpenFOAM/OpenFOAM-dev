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

#include "addToRunTimeSelectionTable.H"

#include "rhoThermo.H"
#include "rhoReactionThermo.H"

#include "rhoCombustionModel.H"

#include "phaseModel.H"
#include "ThermoPhaseModel.H"
#include "IsothermalPhaseModel.H"
#include "AnisothermalPhaseModel.H"
#include "PurePhaseModel.H"
#include "MultiComponentPhaseModel.H"
#include "InertPhaseModel.H"
#include "ReactingPhaseModel.H"
#include "MovingPhaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        MovingPhaseModel
        <
            AnisothermalPhaseModel
            <
                PurePhaseModel
                <
                    InertPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoThermo>
                    >
                >
            >
        >
        purePhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        purePhaseModel,
        phaseSystem,
        purePhaseModel
    );

    typedef
        MovingPhaseModel
        <
            IsothermalPhaseModel
            <
                PurePhaseModel
                <
                    InertPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoThermo>
                    >
                >
            >
        >
        pureIsothermalPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        pureIsothermalPhaseModel,
        phaseSystem,
        pureIsothermalPhaseModel
    );

    typedef
        MovingPhaseModel
        <
            AnisothermalPhaseModel
            <
                MultiComponentPhaseModel
                <
                    InertPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoReactionThermo>
                    >
                >
            >
        >
        multiComponentPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        multiComponentPhaseModel,
        phaseSystem,
        multiComponentPhaseModel
    );

    typedef
        MovingPhaseModel
        <
            AnisothermalPhaseModel
            <
                MultiComponentPhaseModel
                <
                    ReactingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoReactionThermo>,
                        combustionModels::rhoCombustionModel
                    >
                >
            >
        >
        reactingPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        reactingPhaseModel,
        phaseSystem,
        reactingPhaseModel
    );
}

// ************************************************************************* //
