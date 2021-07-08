/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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

#include "fieldsExpression.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldsExpression, 0);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::fieldsExpression::setResultName
(
    const word& functionName,
    const wordList& defaultFieldNames
)
{
    if (fieldNames_.empty())
    {
        fieldNames_ = defaultFieldNames;
    }

    if (resultName_.empty())
    {
        if (!fieldNames_.empty())
        {
            resultName_ = functionName + '(' + fieldNames_[0];
            for (label i=1; i<fieldNames_.size(); i++)
            {
                resultName_ += ',' + fieldNames_[i];
            }
            resultName_ += ')';
        }
        else
        {
            resultName_ = functionName;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldsExpression::fieldsExpression
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const wordList& fieldNames,
    const word& resultName
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldNames_(fieldNames),
    resultName_(resultName)
{
    read(dict);

    if (fieldNames_.size() < 2)
    {
        FatalIOErrorInFunction(dict)
            << "functionObject::" << type() << " " << name
            << " requires at least 2 fields only "
            << fieldNames_.size() << " provided: " << fieldNames_
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldsExpression::~fieldsExpression()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldsExpression::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (fieldNames_.empty() || dict.found("fields"))
    {
        dict.lookup("fields") >> fieldNames_;
    }

    if (dict.found("result"))
    {
        dict.lookup("result") >> resultName_;
    }

    return true;
}


bool Foam::functionObjects::fieldsExpression::execute()
{
    if (!calc())
    {
        DynamicList<word> notFoundFieldNames;
        forAll(fieldNames_, i)
        {
            bool found = false;

            #define findFieldType(Type, GeoField)                              \
                found =                                                        \
                    found                                                      \
                 || mesh_.foundObject<GeoField<Type>>(fieldNames_[i]);
            FOR_ALL_FIELD_TYPES(findFieldType, VolField);
            FOR_ALL_FIELD_TYPES(findFieldType, SurfaceField);
            #undef findFieldType

            if (!found)
            {
                notFoundFieldNames.append(fieldNames_[i]);
            }
        }

        if (!notFoundFieldNames.empty())
        {
            Warning
                << "functionObjects::" << type() << " " << name()
                << " cannot find fields " << notFoundFieldNames << endl;
        }
        else
        {
            Warning
                << "functionObjects::" << type() << " " << name()
                << " fields are not compatible with the " << type()
                << " function" << endl;
        }

        // Clear the result fields from the objectRegistry if present
        clear();

        return false;
    }
    else
    {
        return true;
    }
}


bool Foam::functionObjects::fieldsExpression::write()
{
    return writeObject(resultName_);
}


bool Foam::functionObjects::fieldsExpression::clear()
{
    return clearObject(resultName_);
}


// ************************************************************************* //
