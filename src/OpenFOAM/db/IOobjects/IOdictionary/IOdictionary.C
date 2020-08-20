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

#include "IOdictionary.H"
#include "objectRegistry.H"
#include "Pstream.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IOdictionary, 0);

    bool IOdictionary::writeDictionaries
    (
        debug::infoSwitch("writeDictionaries", 0)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOdictionary::IOdictionary
(
    const IOobject& io,
    const word& wantedType
)
:
    regIOobject(io)
{
    dictionary::name() = IOobject::objectPath();

    // Reading performed by derived type
}


Foam::IOdictionary::IOdictionary(const IOobject& io)
:
    regIOobject(io)
{
    dictionary::name() = IOobject::objectPath();

    readHeaderOk(IOstream::ASCII, typeName);

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


Foam::IOdictionary::IOdictionary
(
    const IOobject& io,
    const dictionary& dict
)
:
    regIOobject(io)
{
    dictionary::name() = IOobject::objectPath();

    if (!readHeaderOk(IOstream::ASCII, typeName))
    {
        dictionary::operator=(dict);
    }

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


Foam::IOdictionary::IOdictionary
(
    const IOobject& io,
    Istream& is
)
:
    regIOobject(io)
{
    dictionary::name() = IOobject::objectPath();

    // Note that we do construct the dictionary null and read in
    // afterwards
    // so that if there is some fancy massaging due to a
    // functionEntry in
    // the dictionary at least the type information is already complete.
    is  >> *this;

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


Foam::IOdictionary::IOdictionary(const IOdictionary& dict)
:
    regIOobject(dict),
    dictionary(dict)
{}


Foam::IOdictionary::IOdictionary(IOdictionary&& dict)
:
    regIOobject(move(dict)),
    dictionary(move(dict))
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::IOdictionary::~IOdictionary()
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::IOdictionary::operator=(const IOdictionary& rhs)
{
    dictionary::operator=(rhs);
}


void Foam::IOdictionary::operator=(IOdictionary&& rhs)
{
    dictionary::operator=(move(rhs));
}


// ************************************************************************* //
