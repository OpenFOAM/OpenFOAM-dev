/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeObjects::writeObjects
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    obr_
    (
        runTime.lookupObject<objectRegistry>
        (
            dict.lookupOrDefault("region", polyMesh::defaultRegion)
        )
    ),
    exclusiveWriting_(false),
    objectNames_()
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
        objectNames_.setSize(1);
        dict.lookup("field") >> objectNames_[0];
    }
    else if (dict.found("fields"))
    {
        dict.lookup("fields") >> objectNames_;
    }
    else
    {
        dict.lookup("objects") >> objectNames_;
    }

    dict.readIfPresent("exclusiveWriting", exclusiveWriting_);

    return true;
}


bool Foam::functionObjects::writeObjects::execute()
{
    return true;
}


bool Foam::functionObjects::writeObjects::write()
{
    Info<< type() << " " << name() << " write:" << nl;

    if (!obr_.time().writeTime())
    {
        obr_.time().writeTimeDict();
    }

    DynamicList<word> allNames(obr_.toc().size());
    forAll(objectNames_, i)
    {
        wordList names(obr_.names<regIOobject>(objectNames_[i]));

        if (names.size())
        {
            allNames.append(names);
        }
        else
        {
            WarningInFunction
                << "Object " << objectNames_[i] << " not found in "
                << "database. Available objects:" << nl << obr_.sortedToc()
                << endl;
        }
    }

    forAll(allNames, i)
    {
        regIOobject& obj = const_cast<regIOobject&>
        (
            obr_.lookupObject<regIOobject>(allNames[i])
        );

        if (exclusiveWriting_)
        {
            // Switch off automatic writing to prevent double write
            obj.writeOpt() = IOobject::NO_WRITE;
        }

        Info<< "    writing object " << obj.name() << endl;

        obj.write();
    }

    return true;
}


// ************************************************************************* //
