/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#define makeMomentumTransportModelTypes(Alpha, Rho, BaseModel)                 \
    namespace Foam                                                             \
    {                                                                          \
        typedef laminarModel<BaseModel> laminar##BaseModel;                    \
        typedef RASModel<BaseModel> RAS##BaseModel;                            \
        typedef LESModel<BaseModel> LES##BaseModel;                            \
    }


#define makeBaseMomentumTransportModel(Alpha, Rho, BaseModel)                  \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef BaseModel Transport##momentumTransportModel;                   \
                                                                               \
        defineRunTimeSelectionTable                                            \
        (                                                                      \
            Transport##momentumTransportModel,                                 \
            dictionary                                                         \
        );                                                                     \
                                                                               \
                                                                               \
        defineNamedTemplateTypeNameAndDebug(laminar##BaseModel, 0);            \
                                                                               \
        defineTemplateRunTimeSelectionTable(laminar##BaseModel, dictionary);   \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            Transport##momentumTransportModel,                                 \
            laminar##BaseModel,                                                \
            dictionary                                                         \
        );                                                                     \
                                                                               \
                                                                               \
        defineNamedTemplateTypeNameAndDebug(RAS##BaseModel, 0);                \
                                                                               \
        defineTemplateRunTimeSelectionTable(RAS##BaseModel, dictionary);       \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            Transport##momentumTransportModel,                                 \
            RAS##BaseModel,                                                    \
            dictionary                                                         \
        );                                                                     \
                                                                               \
                                                                               \
        defineNamedTemplateTypeNameAndDebug(LES##BaseModel, 0);                \
                                                                               \
        defineTemplateRunTimeSelectionTable(LES##BaseModel, dictionary);       \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            Transport##momentumTransportModel,                                 \
            LES##BaseModel,                                                    \
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
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::SType##Models::Type<Foam::BaseModel>, 0);                       \
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
