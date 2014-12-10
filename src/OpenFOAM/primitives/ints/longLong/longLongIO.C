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

Description
    Reads a long long from an input stream, for a given version
    number and File format. If an ascii File is being read, then the line
    numbers are counted and an erroneous read ised.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "longLong.H"
#include "IOstreams.H"

#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::name(long long val)
{
    std::ostringstream buf;
    buf << val;
    return buf.str();
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, long long& l)
{
    l = readLongLong(is);

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, long long&)");

    return is;
}


long long Foam::readLongLong(Istream& is)
{
    register long long result = 0;

    char c = 0;

    static const label zeroOffset = int('0');

    // Get next non-whitespace character
    while (is.read(c) && isspace(c))
    {}

    do
    {
        if (isspace(c) || c == 0) break;

        if (!isdigit(c))
        {
            FatalIOErrorIn("readLongLong(ISstream& is)", is)
                << "Illegal digit: \"" << c << "\""
                << exit(FatalIOError);
        }

        result *= 10 + int(c) - zeroOffset;
    } while (is.read(c));

    return result;
}


bool Foam::readLongLong(const char* buf, long long& s)
{
    char *endptr = NULL;
    s = strtoll(buf, &endptr, 10);
    return (*endptr == 0);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const long long l)
{
    long long val = l;

    long long mask = 1000000000000000000LL;

    bool printZeroes = false;

    while (mask > 0)
    {
        int d = int(val/mask);

        if (d == 0)
        {
            if (printZeroes)
            {
                os.write('0');
            }
        }
        else
        {
            printZeroes = true;
            os.write(char(d+'0'));
        }

        val = val % mask;
        mask /= 10;
    }
    os.check("Ostream& operator<<(Ostream&, const long long)");
    return os;
}


// ************************************************************************* //
