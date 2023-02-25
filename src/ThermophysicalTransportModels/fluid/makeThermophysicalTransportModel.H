/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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

#define makeThermophysicalTransportModelBaseType(                              \
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
        typedef BaseModel<MomentumTransportModel, ThermoModel>                 \
            ThermoModel##MomentumTransportModel##BaseModel;                    \
    }


#define makeThermophysicalTransportModelType(                                  \
    BaseModel, MomentumTransportModel, ThermoModel, SType)                     \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef SType##ThermophysicalTransportModel                            \
        <                                                                      \
             ThermoModel##MomentumTransportModel##BaseModel                    \
        > SType##ThermoModel##MomentumTransportModel##BaseModel;               \
    }


#define makeThermophysicalTransportModelTypes(                                 \
    BaseModel, MomentumTransportModel, ThermoModel)                            \
                                                                               \
    makeThermophysicalTransportModelBaseType                                   \
    (BaseModel, MomentumTransportModel, ThermoModel)                           \
                                                                               \
    makeThermophysicalTransportModelType                                       \
    (BaseModel, MomentumTransportModel, ThermoModel, laminar)                  \
                                                                               \
    makeThermophysicalTransportModelType                                       \
    (BaseModel, MomentumTransportModel, ThermoModel, RAS)                      \
                                                                               \
    makeThermophysicalTransportModelType                                       \
    (BaseModel, MomentumTransportModel, ThermoModel, LES)


#define makeThermophysicalTransportModelTable(                                 \
    BaseModel, MomentumTransportModel, ThermoModel)                            \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            ThermophysicalTransportModel##MomentumTransportModel,              \
            dictionary                                                         \
        );                                                                     \
    }


#define makeThermophysicalTransportModelTableType(                             \
    BaseModel, MomentumTransportModel, ThermoModel, SType)                     \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineNamedTemplateTypeNameAndDebug                                    \
        (                                                                      \
            SType##ThermoModel##MomentumTransportModel##BaseModel,             \
            0                                                                  \
        );                                                                     \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            SType##ThermoModel##MomentumTransportModel##BaseModel,             \
            dictionary                                                         \
        );                                                                     \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            ThermophysicalTransportModel##MomentumTransportModel,              \
            SType##ThermoModel##MomentumTransportModel##BaseModel,             \
            dictionary                                                         \
        );                                                                     \
    }


#define makeThermophysicalTransportModelTables(                                \
    BaseModel, MomentumTransportModel, ThermoModel)                            \
                                                                               \
    makeThermophysicalTransportModelTable                                      \
    (BaseModel, MomentumTransportModel, ThermoModel)                           \
                                                                               \
    makeThermophysicalTransportModelTableType                                  \
    (BaseModel, MomentumTransportModel, ThermoModel, laminar)                  \
                                                                               \
    makeThermophysicalTransportModelTableType                                  \
    (BaseModel, MomentumTransportModel, ThermoModel, RAS)                      \
                                                                               \
    makeThermophysicalTransportModelTableType                                  \
    (BaseModel, MomentumTransportModel, ThermoModel, LES)


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
