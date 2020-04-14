/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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

#define makeMomentumTransportModelTypes(                                       \
    Alpha, Rho, baseModel, BaseModel, Transport)                               \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef BaseModel<Transport>                                           \
            Transport##BaseModel;                                              \
        typedef laminarModel<Transport##BaseModel>                             \
            laminar##Transport##BaseModel;                                     \
        typedef RASModel<Transport##BaseModel>                                 \
            RAS##Transport##BaseModel;                                         \
        typedef LESModel<Transport##BaseModel>                                 \
            LES##Transport##BaseModel;                                         \
    }


#define makeBaseMomentumTransportModel(                                        \
    Alpha, Rho, baseModel, BaseModel, Transport)                               \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef MomentumTransportModel                                         \
        <                                                                      \
            Alpha,                                                             \
            Rho,                                                               \
            baseModel,                                                         \
            Transport                                                          \
        > Transport##baseModel;                                                \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            Transport##baseModel,                                              \
            dictionary                                                         \
        );                                                                     \
                                                                               \
        typedef BaseModel<Transport>                                           \
            Transport##BaseModel;                                              \
                                                                               \
                                                                               \
        typedef laminarModel<Transport##BaseModel>                             \
            laminar##Transport##BaseModel;                                     \
                                                                               \
        defineNamedTemplateTypeNameAndDebug(laminar##Transport##BaseModel, 0); \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (laminar##Transport##BaseModel, dictionary);                           \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            Transport##baseModel,                                              \
            laminar##Transport##BaseModel,                                     \
            dictionary                                                         \
        );                                                                     \
                                                                               \
                                                                               \
        typedef RASModel<Transport##BaseModel>                                 \
            RAS##Transport##BaseModel;                                         \
                                                                               \
        defineNamedTemplateTypeNameAndDebug(RAS##Transport##BaseModel, 0);     \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (RAS##Transport##BaseModel, dictionary);                               \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            Transport##baseModel,                                              \
            RAS##Transport##BaseModel,                                         \
            dictionary                                                         \
        );                                                                     \
                                                                               \
                                                                               \
        typedef LESModel<Transport##BaseModel>                                 \
            LES##Transport##BaseModel;                                         \
                                                                               \
        defineNamedTemplateTypeNameAndDebug(LES##Transport##BaseModel, 0);     \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (LES##Transport##BaseModel, dictionary);                               \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            Transport##baseModel,                                              \
            LES##Transport##BaseModel,                                         \
            dictionary                                                         \
        );                                                                     \
    }


#define makeTemplatedLaminarModel(BaseModel, SType, Type)                      \
    typedef Foam::SType##Models::Type<Foam::BaseModel>                         \
        Type##SType##BaseModel;                                                \
    defineNamedTemplateTypeNameAndDebug(Type##SType##BaseModel, 0);            \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace SType##Models                                                \
        {                                                                      \
            typedef Type<BaseModel> Type##SType##BaseModel;                    \
                                                                               \
            addToRunTimeSelectionTable                                         \
            (                                                                  \
                SType##BaseModel,                                              \
                Type##SType##BaseModel,                                        \
                dictionary                                                     \
            );                                                                 \
        }                                                                      \
    }


#define makeTemplatedMomentumTransportModel(BaseModel, SType, Type)            \
    typedef Foam::SType##Models::Type<Foam::BaseModel>                         \
        Type##SType##BaseModel;                                                \
    defineNamedTemplateTypeNameAndDebug(Type##SType##BaseModel, 0);            \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace SType##Models                                                \
        {                                                                      \
            typedef Type<BaseModel> Type##SType##BaseModel;                    \
                                                                               \
            addToRunTimeSelectionTable                                         \
            (                                                                  \
                SType##BaseModel,                                              \
                Type##SType##BaseModel,                                        \
                dictionary                                                     \
            );                                                                 \
        }                                                                      \
    }


#define makeMomentumTransportModel(BaseModel, SType, Type)                     \
    namespace Foam                                                             \
    {                                                                          \
        namespace SType##Models                                                \
        {                                                                      \
            defineTypeNameAndDebug(Type, 0);                                   \
                                                                               \
            addToRunTimeSelectionTable                                         \
            (                                                                  \
                SType##BaseModel,                                              \
                Type,                                                          \
                dictionary                                                     \
            );                                                                 \
        }                                                                      \
    }


// ************************************************************************* //
