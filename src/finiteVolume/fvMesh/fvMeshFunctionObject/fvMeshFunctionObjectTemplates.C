/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class FieldType>
bool Foam::functionObjects::fvMeshFunctionObject::foundField
(
    const word& fieldName
) const
{
    return mesh_.foundObject<FieldType>(fieldName);
}


template<class FieldType>
const FieldType& Foam::functionObjects::fvMeshFunctionObject::lookupField
(
    const word& fieldName
) const
{
    return mesh_.lookupObject<FieldType>(fieldName);
}


template<class FieldType>
bool Foam::functionObjects::fvMeshFunctionObject::store
(
    word& fieldName,
    const tmp<FieldType>& tfield,
    bool cacheable
)
{
    if (cacheable && fieldName == tfield().name())
    {
        WarningInFunction
            << "Cannot store cache-able field with the named used in the cache."
            << nl
            << "    Either choose a different name or cache the field"
            << "    and use the 'writeRegisteredObject' functionObject."
            << endl;

        return false;
    }

    if
    (
        fieldName.size()
     && mesh_.foundObject<FieldType>(fieldName)
    )
    {
        const_cast<FieldType&>
        (
            mesh_.lookupObject<FieldType>(fieldName)
        ) = tfield;
    }
    else
    {
        if (fieldName.size() && fieldName != tfield().name())
        {
            tfield.ref().rename(fieldName);
        }
        else
        {
            fieldName = tfield().name();
        }

        mesh_.objectRegistry::store(tfield.ptr());
    }

    return true;
}


// ************************************************************************* //
