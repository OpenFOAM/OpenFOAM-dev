/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeThermophysicalTransportModelTypes(                                 \
    BaseModel, MomentumTransportModel, ThermoModel)                            \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef BaseModel<MomentumTransportModel, ThermoModel>                 \
            ThermoModel##MomentumTransportModel##BaseModel;                    \
                                                                               \
        typedef laminarThermophysicalTransportModel                            \
        <                                                                      \
             ThermoModel##MomentumTransportModel##BaseModel                    \
        > laminar##ThermoModel##MomentumTransportModel##BaseModel;             \
                                                                               \
        typedef RASThermophysicalTransportModel                                \
        <                                                                      \
            ThermoModel##MomentumTransportModel##BaseModel                     \
        > RAS##ThermoModel##MomentumTransportModel##BaseModel;                 \
                                                                               \
        typedef LESThermophysicalTransportModel                                \
        <                                                                      \
            ThermoModel##MomentumTransportModel##BaseModel                     \
        > LES##ThermoModel##MomentumTransportModel##BaseModel;                 \
    }


#define makeThermophysicalTransportModels(                                     \
    BaseModel, MomentumTransportModel, ThermoModel)                            \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef ThermophysicalTransportModel                                   \
        <                                                                      \
            MomentumTransportModel,                                            \
            ThermoModel                                                        \
        > ThermophysicalTransportModel##MomentumTransportModel;                \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            ThermophysicalTransportModel##MomentumTransportModel,              \
            dictionary                                                         \
        );                                                                     \
                                                                               \
                                                                               \
        defineNamedTemplateTypeNameAndDebug                                    \
        (                                                                      \
            laminar##ThermoModel##MomentumTransportModel##BaseModel,           \
            0                                                                  \
        );                                                                     \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            laminar##ThermoModel##MomentumTransportModel##BaseModel,           \
            dictionary                                                         \
        );                                                                     \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            ThermophysicalTransportModel##MomentumTransportModel,              \
            laminar##ThermoModel##MomentumTransportModel##BaseModel,           \
            dictionary                                                         \
        );                                                                     \
                                                                               \
                                                                               \
        defineNamedTemplateTypeNameAndDebug                                    \
        (                                                                      \
            RAS##ThermoModel##MomentumTransportModel##BaseModel,               \
            0                                                                  \
        );                                                                     \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            RAS##ThermoModel##MomentumTransportModel##BaseModel,               \
            dictionary                                                         \
        );                                                                     \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            ThermophysicalTransportModel##MomentumTransportModel,              \
            RAS##ThermoModel##MomentumTransportModel##BaseModel,               \
            dictionary                                                         \
        );                                                                     \
                                                                               \
                                                                               \
        defineNamedTemplateTypeNameAndDebug                                    \
        (                                                                      \
            LES##ThermoModel##MomentumTransportModel##BaseModel,               \
            0                                                                  \
        );                                                                     \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            LES##ThermoModel##MomentumTransportModel##BaseModel,               \
            dictionary                                                         \
        );                                                                     \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            ThermophysicalTransportModel##MomentumTransportModel,              \
            LES##ThermoModel##MomentumTransportModel##BaseModel,               \
            dictionary                                                         \
        );                                                                     \
    }


#define makeThermophysicalTransportModel(                                      \
    BaseModel, MomentumTransportModel, ThermoModel, SType, Type)               \
                                                                               \
    typedef Foam::SType##ThermophysicalTransportModels::Type                   \
    <                                                                          \
        Foam::BaseModel<Foam::MomentumTransportModel, Foam::ThermoModel>       \
    > Type##SType##ThermoModel##ThermoModel##MomentumTransportModel##BaseModel;\
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Type##SType##ThermoModel##ThermoModel                                  \
        ##MomentumTransportModel##BaseModel,                                   \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace SType##ThermophysicalTransportModels                         \
        {                                                                      \
            typedef Type<BaseModel<MomentumTransportModel, ThermoModel>>       \
                Type##SType##ThermoModel##ThermoModel                          \
                ##MomentumTransportModel##BaseModel;                           \
                                                                               \
            addToRunTimeSelectionTable                                         \
            (                                                                  \
                SType##ThermoModel##MomentumTransportModel##BaseModel,         \
                Type##SType##ThermoModel##ThermoModel                          \
                ##MomentumTransportModel##BaseModel,                           \
                dictionary                                                     \
            );                                                                 \
        }                                                                      \
    }


#define makeTurbulenceThermophysicalTransportModel(                            \
    BaseModel, MomentumTransportModel, ThermoModel, SType, Type)               \
                                                                               \
    typedef Foam::turbulenceThermophysicalTransportModels::Type                \
    <                                                                          \
        Foam::SType##ThermophysicalTransportModel                              \
        <                                                                      \
            Foam::BaseModel<Foam::MomentumTransportModel, Foam::ThermoModel>   \
        >                                                                      \
    > Type##SType##ThermoModel##ThermoModel##MomentumTransportModel##BaseModel;\
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Type##SType##ThermoModel##ThermoModel                                  \
        ##MomentumTransportModel##BaseModel,                                   \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace turbulenceThermophysicalTransportModels                      \
        {                                                                      \
            typedef Type                                                       \
            <                                                                  \
                SType##ThermophysicalTransportModel                            \
                <                                                              \
                    BaseModel<MomentumTransportModel, ThermoModel>             \
                >                                                              \
            > Type##SType##ThermoModel##MomentumTransportModel##BaseModel;     \
                                                                               \
            addToRunTimeSelectionTable                                         \
            (                                                                  \
                SType##ThermoModel##MomentumTransportModel##BaseModel,         \
                Type##SType##ThermoModel##MomentumTransportModel##BaseModel,   \
                dictionary                                                     \
            );                                                                 \
        }                                                                      \
    }


// ************************************************************************* //
