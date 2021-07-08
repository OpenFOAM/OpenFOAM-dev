/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "uniform.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(uniform, 0);
    addToRunTimeSelectionTable(functionObject, uniform, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::uniform::uniform
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldType_(word::null),
    name_(word::null),
    dimensions_(dimless)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::uniform::~uniform()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::uniform::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    fieldType_ = dict.lookup<word>("fieldType");
    name_ = dict.lookup<word>("name");
    dimensions_.reset(dict.lookup<dimensionSet>("dimensions"));

    bool ok = false;
    #define readValueType(Type, GeoField)                                      \
        if (GeoField<Type>::typeName == fieldType_)                            \
        {                                                                      \
            ok = true;                                                         \
            Type##Value_ = dict.lookup<Type>("value");                         \
        }
    FOR_ALL_FIELD_TYPES(readValueType, VolField);
    FOR_ALL_FIELD_TYPES(readValueType, SurfaceField);
    #undef readValueType

    if (!ok)
    {
        FatalErrorInFunction
            << "Field type " << fieldType_ << " not recognised" << endl << endl;

        DynamicList<word> fieldTypes;
        #define getFieldType(Type, GeoField)                                   \
            fieldTypes.append(GeoField<Type>::typeName);
        FOR_ALL_FIELD_TYPES(getFieldType, VolField);
        FOR_ALL_FIELD_TYPES(getFieldType, SurfaceField);
        #undef getFieldType

        FatalErrorInFunction
            << "Available field types are : " << endl
            << fieldTypes
            << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::uniform::execute()
{
    #define calcType(Type, GeoField)                                           \
        if (GeoField<Type>::typeName == fieldType_)                            \
        {                                                                      \
            store                                                              \
            (                                                                  \
                name_,                                                         \
                GeoField<Type>::New                                            \
                (                                                              \
                    name_,                                                     \
                    mesh_,                                                     \
                    dimensioned<Type>(dimensions_, Type##Value_)               \
                )                                                              \
            );                                                                 \
        }
    FOR_ALL_FIELD_TYPES(calcType, VolField);
    FOR_ALL_FIELD_TYPES(calcType, SurfaceField);
    #undef calcType

    return true;
}


bool Foam::functionObjects::uniform::write()
{
    return writeObject(name_);
}


bool Foam::functionObjects::uniform::clear()
{
    return clearObject(name_);
}


// ************************************************************************* //
