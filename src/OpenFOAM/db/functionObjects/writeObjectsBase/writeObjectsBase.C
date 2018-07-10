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

#include "writeObjectsBase.H"
#include "Time.H"
#include "dictionary.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::writeObjectsBase::resetWriteObjectName
(
    const wordRe& name
)
{
    writeObjectNames_.clear();
    writeObjectNames_.append(name);
}


void Foam::functionObjects::writeObjectsBase::resetWriteObjectNames
(
    const wordReList& names
)
{
    writeObjectNames_.clear();
    writeObjectNames_.append(names);
}


Foam::wordList Foam::functionObjects::writeObjectsBase::objectNames()
{
    DynamicList<word> allNames(writeObr_.toc().size());
    forAll(writeObjectNames_, i)
    {
        wordList names(writeObr_.names<regIOobject>(writeObjectNames_[i]));

        if (names.size())
        {
            allNames.append(names);
        }
        else
        {
            WarningInFunction
                << "Object " << writeObjectNames_[i] << " not found in "
                << "database. Available objects:" << nl << writeObr_.sortedToc()
                << endl;
        }
    }

    return allNames;
}


void Foam::functionObjects::writeObjectsBase::writeObject
(
    const regIOobject& obj
)
{
    if (log_) Info << "    writing object " << obj.name() << endl;

    obj.write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeObjectsBase::writeObjectsBase
(
    const objectRegistry& obr,
    const Switch& log
)
:
    writeObr_(obr),
    log_(log)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::writeObjectsBase::~writeObjectsBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordReList&
Foam::functionObjects::writeObjectsBase::writeObjectNames() const
{
    return writeObjectNames_;
}


bool Foam::functionObjects::writeObjectsBase::read(const dictionary& dict)
{
    dict.lookup("objects") >> writeObjectNames_;

    return true;
}


bool Foam::functionObjects::writeObjectsBase::write()
{
    wordList names(objectNames());

    if(!names.empty())
    {
        if (!writeObr_.time().writeTime())
        {
            writeObr_.time().writeTimeDict();
        }

        forAll(names, i)
        {
            const regIOobject& obj =
                writeObr_.lookupObject<regIOobject>(names[i]);

            writeObject(obj);
        }
    }

    return true;
}


// ************************************************************************* //
