/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "regionFunctionObject.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(regionFunctionObject, 0);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::functionObjects::regionFunctionObject::writeObject
(
    const word& fieldName
)
{
    if (obr_.foundObject<regIOobject>(fieldName))
    {
        const regIOobject& field = obr_.lookupObject<regIOobject>(fieldName);

        Log << "    functionObjects::" << type() << " " << name()
            << " writing field: " << field.name() << endl;

        field.write();

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::functionObjects::regionFunctionObject::clearObject
(
    const word& fieldName
)
{
    if (foundObject<regIOobject>(fieldName))
    {
        regIOobject& resultObject = lookupObjectRef<regIOobject>(fieldName);

        if (resultObject.ownedByRegistry())
        {
            return resultObject.checkOut();
        }
        else
        {
            return false;
        }
    }
    else
    {
        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::regionFunctionObject::regionFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime),
    obr_
    (
        runTime.lookupObject<objectRegistry>
        (
            dict.lookupOrDefault("region", polyMesh::defaultRegion)
        )
    )
{}


Foam::functionObjects::regionFunctionObject::regionFunctionObject
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    functionObject(name),
    time_(obr.time()),
    obr_(obr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::regionFunctionObject::~regionFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::regionFunctionObject::read(const dictionary& dict)
{
    return functionObject::read(dict);
}


// ************************************************************************* //
