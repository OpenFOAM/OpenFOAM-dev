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

#include "writeLocalObjects.H"
#include "stringListOps.H"
#include "dictionary.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::writeLocalObjects::resetLocalObjectName
(
    const word& name
)
{
    localObjectNames_.clear();
    localObjectNames_.append(name);
}


void Foam::functionObjects::writeLocalObjects::resetLocalObjectNames
(
    const wordList& names
)
{
    localObjectNames_.clear();
    localObjectNames_.append(names);
}


Foam::wordList Foam::functionObjects::writeLocalObjects::objectNames()
{
    wordList names
    (
        subsetStrings(wordReListMatcher(writeObjectNames_), localObjectNames_)
    );

    return names;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeLocalObjects::writeLocalObjects
(
    const objectRegistry& obr,
    const Switch& log
)
:
    writeObjectsBase(obr, log),
    localObjectNames_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::writeLocalObjects::~writeLocalObjects()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList&
Foam::functionObjects::writeLocalObjects::localObjectNames() const
{
    return localObjectNames_;
}


bool Foam::functionObjects::writeLocalObjects::read(const dictionary& dict)
{
    if (dict.found("objects"))
    {
        writeObjectsBase::read(dict);
    }
    else
    {
        resetWriteObjectName(wordRe(".*", wordRe::compOption::detect));
    }

    return true;
}


// ************************************************************************* //
