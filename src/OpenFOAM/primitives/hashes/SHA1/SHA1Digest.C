/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "SHA1Digest.H"
#include "IOstreams.H"

#include <cstring>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::SHA1Digest Foam::SHA1Digest::null;

//! \cond fileScope
static const char hexChars[] = "0123456789abcdef";
//! \endcond


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

unsigned char Foam::SHA1Digest::readHexDigit(Istream& is)
{
    // Takes into account that 'a' (or 'A') is 10
    static const int alphaOffset = toupper('A') - 10;
    // Takes into account that '0' is 0
    static const int zeroOffset = int('0');


    // silently ignore leading or intermediate '_'
    char c = 0;
    do
    {
        is.read(c);
    }
    while (c == '_');

    if (!isxdigit(c))
    {
        FatalIOErrorIn("SHA1Digest::readHexDigit(Istream&)", is)
            << "Illegal hex digit: '" << c << "'"
            << exit(FatalIOError);
    }

    if (isdigit(c))
    {
        return int(c) - zeroOffset;
    }
    else
    {
        return toupper(c) - alphaOffset;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SHA1Digest::SHA1Digest()
{
    clear();
}


Foam::SHA1Digest::SHA1Digest(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::SHA1Digest::clear()
{
    memset(v_, 0, length);
}


bool Foam::SHA1Digest::empty() const
{
    for (unsigned i = 0; i < length; ++i)
    {
        if (v_[i])
        {
            return false;
        }
    }

    return true;
}


std::string Foam::SHA1Digest::str(const bool prefixed) const
{
    std::string buf;
    unsigned nChar = 0;

    if (prefixed)
    {
        buf.resize(1 + length*2);
        buf[nChar++] = '_';
    }
    else
    {
        buf.resize(length*2);
    }

    for (unsigned i = 0; i < length; ++i)
    {
        buf[nChar++] = hexChars[((v_[i] >> 4) & 0xF)];
        buf[nChar++] = hexChars[(v_[i] & 0xF)];
    }

    return buf;
}


Foam::Ostream& Foam::SHA1Digest::write(Ostream& os, const bool prefixed) const
{
    if (prefixed)
    {
        os.write('_');
    }

    for (unsigned i = 0; i < length; ++i)
    {
        os.write(hexChars[((v_[i] >> 4) & 0xF)]);
        os.write(hexChars[(v_[i] & 0xF)]);
    }

    os.check("SHA1Digest::write(Ostream&, const bool)");
    return os;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

bool Foam::SHA1Digest::operator==(const SHA1Digest& rhs) const
{
    for (unsigned i = 0; i < length; ++i)
    {
        if (v_[i] != rhs.v_[i])
        {
            return false;
        }
    }

    return true;
}


bool Foam::SHA1Digest::operator==(const std::string& hexdigits) const
{
    // null or empty string is not an error - interpret as '0000..'
    if (hexdigits.empty())
    {
        return empty();
    }

    // skip possible '_' prefix
    unsigned charI = 0;
    if (hexdigits[0] == '_')
    {
        ++charI;
    }

    // incorrect length - can never match
    if (hexdigits.size() != charI + length*2)
    {
        return false;
    }

    for (unsigned i = 0; i < length; ++i)
    {
        const char c1 = hexChars[((v_[i] >> 4) & 0xF)];
        const char c2 = hexChars[(v_[i] & 0xF)];

        if (c1 != hexdigits[charI++]) return false;
        if (c2 != hexdigits[charI++]) return false;
    }

    return true;
}


bool Foam::SHA1Digest::operator==(const char* hexdigits) const
{
    // null or empty string is not an error - interpret as '0000..'
    if (!hexdigits || !*hexdigits)
    {
        return empty();
    }

    // skip possible '_' prefix
    unsigned charI = 0;
    if (hexdigits[0] == '_')
    {
        ++charI;
    }

    // incorrect length - can never match
    if (strlen(hexdigits) != charI + length*2)
    {
        return false;
    }

    for (unsigned i = 0; i < length; ++i)
    {
        const char c1 = hexChars[((v_[i] >> 4) & 0xF)];
        const char c2 = hexChars[(v_[i] & 0xF)];

        if (c1 != hexdigits[charI++]) return false;
        if (c2 != hexdigits[charI++]) return false;
    }

    return true;
}


bool Foam::SHA1Digest::operator!=(const SHA1Digest& rhs) const
{
    return !operator==(rhs);
}


bool Foam::SHA1Digest::operator!=(const std::string& rhs) const
{
    return !operator==(rhs);
}


bool Foam::SHA1Digest::operator!=(const char* rhs) const
{
    return !operator==(rhs);
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, SHA1Digest& dig)
{
    unsigned char *v = dig.v_;

    for (unsigned i = 0; i < dig.length; ++i)
    {
        unsigned char c1 = SHA1Digest::readHexDigit(is);
        unsigned char c2 = SHA1Digest::readHexDigit(is);

        v[i] = (c1 << 4) + c2;
    }

    is.check("Istream& operator>>(Istream&, SHA1Digest&)");
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const SHA1Digest& dig)
{
    return dig.write(os);
}


// ************************************************************************* //
