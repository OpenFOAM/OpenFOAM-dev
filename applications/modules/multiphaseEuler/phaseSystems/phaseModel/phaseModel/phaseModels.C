/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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
#include "rhoMulticomponentThermo.H"

#include "combustionModel.H"

#include "phaseModel.H"
#include "ThermoPhaseModel.H"
#include "IsothermalPhaseModel.H"
#include "AnisothermalPhaseModel.H"
#include "PurePhaseModel.H"
#include "MulticomponentPhaseModel.H"
#include "InertPhaseModel.H"
#include "ReactingPhaseModel.H"
#include "MovingPhaseModel.H"
#include "StationaryPhaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        AnisothermalPhaseModel
        <
            PurePhaseModel
            <
                InertPhaseModel
                <
                    MovingPhaseModel
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
        AnisothermalPhaseModel
        <
            PurePhaseModel
            <
                InertPhaseModel
                <
                    StationaryPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoThermo>
                    >
                >
            >
        >
        pureStationaryPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        pureStationaryPhaseModel,
        phaseSystem,
        pureStationaryPhaseModel
    );

    typedef
        IsothermalPhaseModel
        <
            PurePhaseModel
            <
                InertPhaseModel
                <
                    MovingPhaseModel
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
        IsothermalPhaseModel
        <
            PurePhaseModel
            <
                InertPhaseModel
                <
                    StationaryPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoThermo>
                    >
                >
            >
        >
        pureStationaryIsothermalPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        pureStationaryIsothermalPhaseModel,
        phaseSystem,
        pureStationaryIsothermalPhaseModel
    );

    typedef
        AnisothermalPhaseModel
        <
            MulticomponentPhaseModel
            <
                InertPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoMulticomponentThermo>
                    >
                >
            >
        >
        multicomponentPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        multicomponentPhaseModel,
        phaseSystem,
        multicomponentPhaseModel
    );

    // Also add multicomponentPhaseModel with the name multiComponentPhaseModel
    // for backward-compatibility
    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        multicomponentPhaseModel,
        phaseSystem,
        multiComponentPhaseModel
    );

    typedef
        IsothermalPhaseModel
        <
            MulticomponentPhaseModel
            <
                InertPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoMulticomponentThermo>
                    >
                >
            >
        >
        multicomponentIsothermalPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        multicomponentIsothermalPhaseModel,
        phaseSystem,
        multicomponentIsothermalPhaseModel
    );

    typedef
        AnisothermalPhaseModel
        <
            MulticomponentPhaseModel
            <
                ReactingPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoMulticomponentThermo>
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
