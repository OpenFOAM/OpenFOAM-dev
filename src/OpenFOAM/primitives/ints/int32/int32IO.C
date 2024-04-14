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

#include "int32.H"
#include "IOstreams.H"

#include <inttypes.h>
#include <sstream>
#include <cerrno>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::name(const int32_t val)
{
    std::ostringstream buf;
    buf << val;
    return buf.str();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, int32_t& i)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isInteger32())
    {
        i = t.integer32Token();
    }
    else
    {
        is.setBad();
        FatalIOErrorInFunction(is)
            << "wrong token type - expected int32_t, found " << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, int32_t&)");

    return is;
}


int32_t Foam::readInt32(Istream& is)
{
    int32_t val;
    is >> val;

    return val;
}


bool Foam::read(const char* buf, int32_t& s)
{
    char *endptr = nullptr;
    errno = 0;
    intmax_t l = strtoimax(buf, &endptr, 10);
    s = int32_t(l);
    return
        (*endptr == 0) && (errno == 0)
     && (l >= INT32_MIN) && (l <= INT32_MAX);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const int32_t i)
{
    os.write(i);
    os.check("Ostream& operator<<(Ostream&, const int32_t)");
    return os;
}


#if WM_ARCH_OPTION == 32
Foam::Istream& Foam::operator>>(Istream& is, long& i)
{
    return operator>>(is, reinterpret_cast<int32_t&>(i));
}


Foam::Ostream& Foam::operator<<(Ostream& os, const long i)
{
    os << int32_t(i);
    return os;
}
#endif


// ************************************************************************* //
