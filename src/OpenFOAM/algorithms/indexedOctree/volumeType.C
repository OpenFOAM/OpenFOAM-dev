/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "volumeType.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::volumeType,
        4
    >::names[] =
    {
        "unknown",
        "mixed",
        "inside",
        "outside"
    };
}

const Foam::NamedEnum<Foam::volumeType, 4> Foam::volumeType::names;


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, volumeType& vt)
{
    // Read beginning of volumeType
    is.readBegin("volumeType");

    int type;
    is  >> type;

    vt.t_ = static_cast<volumeType::type>(type);

    // Read end of volumeType
    is.readEnd("volumeType");

    // Check state of Istream
    is.check("operator>>(Istream&, volumeType&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const volumeType& vt)
{
    os  << static_cast<int>(vt.t_);

    return os;
}


// ************************************************************************* //
