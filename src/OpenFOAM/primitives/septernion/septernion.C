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

#include "septernion.H"
#include "IOstreams.H"
#include "OStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::septernion::typeName = "septernion";
const Foam::septernion Foam::septernion::zero
(
    vector(0, 0, 0),
    quaternion(0, vector(0, 0, 0))
);
const Foam::septernion Foam::septernion::I
(
    vector(0, 0, 0),
    quaternion(1, vector(0, 0, 0))
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::septernion::septernion(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::name(const septernion& s)
{
    OStringStream buf;
    buf << '(' << s.t() << ',' << s.r() << ')';
    return buf.str();
}


Foam::septernion Foam::slerp
(
    const septernion& qa,
    const septernion& qb,
    const scalar t
)
{
    return septernion((1.0-t)*qa.t()+t*qb.t(), slerp(qa.r(), qb.r(), t));
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, septernion& s)
{
    // Read beginning of septernion
    is.readBegin("septernion");

    is  >> s.t() >> s.r();

    // Read end of septernion
    is.readEnd("septernion");

    // Check state of Istream
    is.check("operator>>(Istream&, septernion&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const septernion& s)
{
    os  << token::BEGIN_LIST
        << s.t() << token::SPACE << s.r()
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
