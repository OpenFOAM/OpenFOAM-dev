/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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

#include "localIOdictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::localIOdictionary::localIOdictionary(const IOobject& io)
:
    baseIOdictionary(io)
{
    readHeaderOk(IOstream::ASCII, typeName);

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


Foam::localIOdictionary::localIOdictionary
(
    const IOobject& io,
    const word& wantedType
)
:
    baseIOdictionary(io)
{
    readHeaderOk(IOstream::ASCII, wantedType);

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


Foam::localIOdictionary::localIOdictionary
(
    const IOobject& io,
    const dictionary& dict
)
:
    baseIOdictionary(io, dict)
{
    if (!readHeaderOk(IOstream::ASCII, typeName))
    {
        dictionary::operator=(dict);
    }

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


Foam::localIOdictionary::localIOdictionary
(
    const IOobject& io,
    Istream& is
)
:
    baseIOdictionary(io, is)
{
    // Note that we do construct the dictionary null and read in
    // afterwards
    // so that if there is some fancy massaging due to a
    // functionEntry in
    // the dictionary at least the type information is already complete.
    is  >> *this;

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


Foam::localIOdictionary::localIOdictionary
(
    localIOdictionary&& dict
)
:
    baseIOdictionary(move(dict))
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::localIOdictionary::~localIOdictionary()
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::localIOdictionary::operator=(localIOdictionary&& rhs)
{
    baseIOdictionary::operator=(move(rhs));
}


// ************************************************************************* //
