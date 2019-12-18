/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "surfZone.H"
#include "dictionary.H"
#include "word.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(surfZone, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfZone::surfZone()
:
    surfZoneIdentifier(),
    size_(0),
    start_(0)
{}


Foam::surfZone::surfZone
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& geometricType
)
:
    surfZoneIdentifier(name, index, geometricType),
    size_(size),
    start_(start)
{}


Foam::surfZone::surfZone(Istream& is, const label index)
:
    surfZoneIdentifier(),
    size_(0),
    start_(0)
{
    word name(is);
    dictionary dict(is);

    operator=(surfZone(name, dict, index));
}


Foam::surfZone::surfZone
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    surfZoneIdentifier(name, dict, index),
    size_(dict.lookup<label>("nFaces")),
    start_(dict.lookup<label>("startFace"))
{}


Foam::surfZone::surfZone(const surfZone& zone, const label index)
:
    surfZoneIdentifier(zone, index),
    size_(zone.size()),
    start_(zone.start())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfZone::write(Ostream& os) const
{
    writeDict(os);
}


void Foam::surfZone::writeDict(Ostream& os) const
{
    os  << indent << name() << nl
        << indent << token::BEGIN_BLOCK << incrIndent << nl;

    surfZoneIdentifier::write(os);
    writeEntry(os, "nFaces", size());
    writeEntry(os, "startFace", start());

    os  << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::surfZone::operator!=(const surfZone& rhs) const
{
    return !(*this == rhs);
}


bool Foam::surfZone::operator==(const surfZone& rhs) const
{
    return
    (
        size() == rhs.size()
     && start() == rhs.start()
     && geometricType() == rhs.geometricType()
    );
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, surfZone& zone)
{
    zone = surfZone(is, 0);

    is.check("Istream& operator>>(Istream&, surfZone&)");
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const surfZone& zone)
{
    zone.write(os);
    os.check("Ostream& operator<<(Ostream&, const surfZone&");
    return os;
}


// ************************************************************************* //
