/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "writeObjects.H"
#include "Time.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(writeObjects, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeObjects,
        dictionary
    );
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::writeObjects::writeOption,
    3
>::names[] =
{
    "autoWrite",
    "noWrite",
    "anyWrite"
};

const Foam::NamedEnum
<
    Foam::functionObjects::writeObjects::writeOption,
    3
> Foam::functionObjects::writeObjects::writeOptionNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::writeObjects::writeObject
(
    const regIOobject& obj
)
{
    switch (writeOption_)
    {
        case writeOption::AUTO_WRITE:
        {
            if (obj.writeOpt() != IOobject::AUTO_WRITE)
            {
                return;
            }

            break;
        }
        case writeOption::NO_WRITE:
        {
            if (obj.writeOpt() != IOobject::NO_WRITE)
            {
                return;
            }

            break;
        }
        case writeOption::ANY_WRITE:
        {
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown writeOption "
                << writeOptionNames_[writeOption_]
                << ". Valid writeOption types are" << writeOptionNames_
                << exit(FatalError);
        }
    }

    if
    (
        obj.writeOpt() == IOobject::AUTO_WRITE
     && writeObr_.time().writeTime()
    )
    {
        Log << "    automatically written object " << obj.name() << endl;
    }
    else
    {
        if (obj.db().cacheTemporaryObject(obj.name()))
        {
            // If the object is a temporary field expression wrap with tmp<...>
            const word name(obj.name());
            regIOobject& objRef = const_cast<regIOobject&>(obj);
            objRef.IOobject::rename("tmp<" + name + ">");
            writeObjectsBase::writeObject(obj);
            objRef.IOobject::rename(name);
        }
        else
        {
            writeObjectsBase::writeObject(obj);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeObjects::writeObjects
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    writeObjectsBase
    (
        runTime.lookupObject<objectRegistry>
        (
            dict.lookupOrDefault("region", polyMesh::defaultRegion)
        ),
        log
    ),
    writeOption_(writeOption::ANY_WRITE)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::writeObjects::~writeObjects()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeObjects::read(const dictionary& dict)
{
    if (dict.found("field"))
    {
        writeObjectNames_.setSize(1);
        dict.lookup("field") >> writeObjectNames_[0];
    }
    else if (dict.found("fields"))
    {
        dict.lookup("fields") >> writeObjectNames_;
    }
    else
    {
        writeObjectsBase::read(dict);
    }

    if (dict.found("writeOption"))
    {
        writeOption_ = writeOptionNames_.read(dict.lookup("writeOption"));
    }
    else
    {
        writeOption_ = writeOption::ANY_WRITE;
    }

    executeAtStart_ = dict.lookupOrDefault<Switch>("executeAtStart", false);

    return functionObject::read(dict);
}


bool Foam::functionObjects::writeObjects::execute()
{
    return true;
}


bool Foam::functionObjects::writeObjects::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeObjectsBase::write();

    Log << endl;

    return true;
}


// ************************************************************************* //
