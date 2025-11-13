/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "OTstream.H"
#include "int.H"
#include "token.H"

#include <cctype>

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::OTstream::OTstream
(
    const string& name,
    const streamFormat format,
    const versionNumber version,
    const bool global
)
:
    Ostream(format, version, UNCOMPRESSED, global),
    name_(name)
{
    setOpened();
    setGood();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::OTstream::~OTstream()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::OTstream::write(const token& t)
{
    append(t);
    return *this;
}


Foam::Ostream& Foam::OTstream::write(const char c)
{
    switch (c)
    {
        // Punctuation
        case token::END_STATEMENT :
        case token::BEGIN_LIST :
        case token::END_LIST :
        case token::BEGIN_SQR :
        case token::END_SQR :
        case token::BEGIN_BLOCK :
        case token::END_BLOCK :
        case token::COLON :
        case token::COMMA :
        case token::ASSIGN :
        case token::ADD :
        case token::SUBTRACT :
        case token::MULTIPLY :
        case token::DIVIDE :
        {
            append(token::punctuationToken(c));
            return *this;
        }

        default:
        {
            if (!isspace(c))
            {
                append(token::punctuationToken(c));
            }
        }
    }

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const char* str)
{
    const word nonWhiteChars(string::validate<word>(str));

    if (nonWhiteChars.size() == 1)
    {
        append(nonWhiteChars[0]);
    }
    else if (nonWhiteChars.size())
    {
        append(nonWhiteChars);
    }

    return *this;
}


Foam::Ostream& Foam::OTstream::write(const word& str)
{
    append(str);
    return *this;
}


Foam::Ostream& Foam::OTstream::write(const string& str)
{
    append(str);
    return *this;
}


Foam::Ostream& Foam::OTstream::write(const keyType& kt)
{
    append(kt);
    return *this;
}


Foam::Ostream& Foam::OTstream::write(const verbatimString& vs)
{
    append(vs);
    return *this;
}


Foam::Ostream& Foam::OTstream::writeQuoted
(
    const std::string& str,
    const bool quoted
)
{
    if (quoted)
    {
        append(string(str));
    }
    else
    {
        append(word(str));
    }
    return *this;
}


Foam::Ostream& Foam::OTstream::write(const int32_t val)
{
    append(val);
    return *this;
}


Foam::Ostream& Foam::OTstream::write(const int64_t val)
{
    append(val);
    return *this;
}


Foam::Ostream& Foam::OTstream::write(const uint32_t val)
{
    append(val);
    return *this;
}


Foam::Ostream& Foam::OTstream::write(const uint64_t val)
{
    append(val);
    return *this;
}


Foam::Ostream& Foam::OTstream::write(const floatScalar val)
{
    append(val);
    return *this;
}


Foam::Ostream& Foam::OTstream::write(const doubleScalar val)
{
    append(val);
    return *this;
}


Foam::Ostream& Foam::OTstream::write(const longDoubleScalar val)
{
    append(val);
    return *this;
}


Foam::Ostream& Foam::OTstream::write(const char* data, std::streamsize count)
{
    NotImplemented;
    return *this;
}


void Foam::OTstream::print(Ostream& os) const
{
    os  << "OTstream: " << name().c_str() << ' ';
    IOstream::print(os);
}


// ************************************************************************* //
