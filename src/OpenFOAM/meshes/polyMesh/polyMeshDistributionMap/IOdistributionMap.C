/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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

#include "IOdistributionMap.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(IOdistributionMap, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOdistributionMap::IOdistributionMap(const IOobject& io)
:
    regIOobject(io)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningInFunction
            << "IOdistributionMap " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOdistributionMap does not support automatic rereading."
            << endl;
    }

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
}


Foam::IOdistributionMap::IOdistributionMap
(
    const IOobject& io,
    const distributionMap& map
)
:
    regIOobject(io)
{
     // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningInFunction
            << "IOdistributionMap " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOdistributionMap does not support automatic rereading."
            << endl;
    }

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
    else
    {
        distributionMap::operator=(map);
    }
}


Foam::IOdistributionMap::IOdistributionMap
(
    const IOobject& io,
    distributionMap&& map
)
:
    regIOobject(io),
    distributionMap(move(map))
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningInFunction
            << "IOdistributionMap " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOdistributionMap does not support automatic rereading."
            << endl;
    }

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::IOdistributionMap::~IOdistributionMap()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::IOdistributionMap::readData(Istream& is)
{
    return (is >> *this).good();
}


bool Foam::IOdistributionMap::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// ************************************************************************* //
