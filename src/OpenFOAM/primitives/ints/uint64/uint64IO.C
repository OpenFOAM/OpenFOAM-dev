/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2024 OpenFOAM Foundation
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

#include "uint64.H"
#include "IOstreams.H"

#include <inttypes.h>
#include <sstream>
#include <cerrno>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::name(const uint64_t val)
{
    std::ostringstream buf;
    buf << val;
    return buf.str();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, uint64_t& i)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isUnsignedInteger64())
    {
        i = t.unsignedInteger64Token();
    }
    else
    {
        is.setBad();
        FatalIOErrorInFunction(is)
            << "wrong token type - expected uint64_t, found " << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, uint64_t&)");

    return is;
}


uint64_t Foam::readUint64(Istream& is)
{
    uint64_t val;
    is >> val;

    return val;
}


bool Foam::read(const char* buf, uint64_t& s)
{
    char *endptr = nullptr;
    errno = 0;
    uintmax_t l = strtoumax(buf, &endptr, 10);
    s = uint64_t(l);
    return (*endptr == 0) && (errno == 0);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const uint64_t i)
{
    os.write(i);
    os.check("Ostream& operator<<(Ostream&, const uint64_t)");
    return os;
}


// ************************************************************************* //
