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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeThermophysicalTransportModelTypes(                                 \
    BaseModel, MomentumTransportModel)                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef BaseModel<MomentumTransportModel>                              \
            BaseModel##MomentumTransportModel;                                 \
                                                                               \
        typedef laminarThermophysicalTransportModel                            \
        <                                                                      \
             BaseModel##MomentumTransportModel                                 \
        > laminar##BaseModel##MomentumTransportModel;                          \
                                                                               \
        typedef RASThermophysicalTransportModel                                \
        <                                                                      \
            BaseModel##MomentumTransportModel                                  \
        > RAS##BaseModel##MomentumTransportModel;                              \
                                                                               \
        typedef LESThermophysicalTransportModel                                \
        <                                                                      \
            BaseModel##MomentumTransportModel                                  \
        > LES##BaseModel##MomentumTransportModel;                              \
    }


#define makeThermophysicalTransportModel(BaseModel, MomentumTransportModel)    \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef ThermophysicalTransportModel                                   \
        <                                                                      \
            MomentumTransportModel                                             \
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
            laminar##BaseModel##MomentumTransportModel,                        \
            0                                                                  \
        );                                                                     \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            laminar##BaseModel##MomentumTransportModel,                        \
            dictionary                                                         \
        );                                                                     \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            ThermophysicalTransportModel##MomentumTransportModel,              \
            laminar##BaseModel##MomentumTransportModel,                        \
            dictionary                                                         \
        );                                                                     \
                                                                               \
                                                                               \
        defineNamedTemplateTypeNameAndDebug                                    \
        (                                                                      \
            RAS##BaseModel##MomentumTransportModel,                            \
            0                                                                  \
        );                                                                     \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            RAS##BaseModel##MomentumTransportModel,                            \
            dictionary                                                         \
        );                                                                     \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            ThermophysicalTransportModel##MomentumTransportModel,              \
            RAS##BaseModel##MomentumTransportModel,                            \
            dictionary                                                         \
        );                                                                     \
                                                                               \
                                                                               \
        defineNamedTemplateTypeNameAndDebug                                    \
        (                                                                      \
            LES##BaseModel##MomentumTransportModel,\
            0                                                                  \
        );                                                                     \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            LES##BaseModel##MomentumTransportModel,\
            dictionary                                                         \
        );                                                                     \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            ThermophysicalTransportModel##MomentumTransportModel,              \
            LES##BaseModel##MomentumTransportModel,\
            dictionary                                                         \
        );                                                                     \
    }


#define makeTemplatedThermophysicalTransportModel(                             \
    BaseModel, MomentumTransportModel, SType, Type)                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::SType##ThermophysicalTransportModels::Type                       \
        <                                                                      \
            Foam::BaseModel<Foam::MomentumTransportModel>                      \
        >,                                                                     \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace SType##ThermophysicalTransportModels                         \
        {                                                                      \
            typedef Type<BaseModel<MomentumTransportModel>>                    \
                Type##SType##BaseModel##MomentumTransportModel;                \
                                                                               \
            addToRunTimeSelectionTable                                         \
            (                                                                  \
                SType##BaseModel##MomentumTransportModel,                      \
                Type##SType##BaseModel##MomentumTransportModel,                \
                dictionary                                                     \
            );                                                                 \
        }                                                                      \
    }


// ************************************************************************* //
