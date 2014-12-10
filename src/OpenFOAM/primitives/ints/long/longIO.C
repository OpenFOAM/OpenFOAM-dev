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
    Reads a long from an input stream, for a given version
    number and File format. If an ascii File is being read, then the line
    numbers are counted and an erroneous read ised.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "long.H"
#include "IOstreams.H"

#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::name(const long val)
{
    std::ostringstream buf;
    buf << val;
    return buf.str();
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, long& l)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isLabel())
    {
        l = long(t.labelToken());
    }
    else
    {
        is.setBad();
        FatalIOErrorIn("operator>>(Istream&, long&)", is)
            << "wrong token type - expected long, found " << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, long&)");

    return is;
}


long Foam::readLong(Istream& is)
{
    long val;
    is >> val;

    return val;
}

bool Foam::readLong(const char* buf, long& s)
{
    char *endptr = NULL;
    s = strtol(buf, &endptr, 10);
    return (*endptr == 0);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const long l)
{
    os.write(label(l));
    os.check("Ostream& operator<<(Ostream&, const long)");
    return os;
}


// ************************************************************************* //
