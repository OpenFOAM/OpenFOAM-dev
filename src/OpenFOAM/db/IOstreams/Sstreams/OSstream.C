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

#include "error.H"
#include "OSstream.H"
#include "token.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Ostream& Foam::OSstream::write(const char c)
{
    os_ << c;
    if (c == token::NL)
    {
        lineNumber_++;
    }
    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const char* str)
{
    lineNumber_ += string(str).count(token::NL);
    os_ << str;
    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const word& str)
{
    os_ << str;
    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const string& str)
{
    return writeQuoted(str);
}


Foam::Ostream& Foam::OSstream::write(const verbatimString& vs)
{
    os_ << token::HASH << token::BEGIN_BLOCK;
    writeQuoted(vs, false);
    os_ << token::HASH << token::END_BLOCK;
    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::writeQuoted
(
    const std::string& str,
    const bool quoted
)
{
    if (quoted)
    {
        os_ << token::BEGIN_STRING;

        int backslash = 0;
        for
        (
            string::const_iterator iter = str.begin();
            iter != str.end();
            ++iter
        )
        {
            char c = *iter;

            if (c == '\\')
            {
                backslash++;
                // suppress output until we know if other characters follow
                continue;
            }
            else if (c == token::NL)
            {
                lineNumber_++;
                backslash++;    // backslash escape for newline
            }
            else if (c == token::END_STRING)
            {
                backslash++;    // backslash escape for quote
            }

            // output pending backslashes
            while (backslash)
            {
                os_ << '\\';
                backslash--;
            }

            os_ << c;
        }

        // silently drop any trailing backslashes
        // they would otherwise appear like an escaped end-quote
        os_ << token::END_STRING;
    }
    else
    {
        // output unquoted string, only advance line number on newline
        lineNumber_ += string(str).count(token::NL);
        os_ << str;
    }

    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const int32_t val)
{
    os_ << val;
    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const int64_t val)
{
    os_ << val;
    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const floatScalar val)
{
    os_ << val;
    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const doubleScalar val)
{
    os_ << val;
    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const longDoubleScalar val)
{
    os_ << val;
    setState(os_.rdstate());
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const char* buf, std::streamsize count)
{
    if (format() != BINARY)
    {
        FatalIOErrorInFunction(*this)
            << "stream format not binary"
            << abort(FatalIOError);
    }

    os_ << token::BEGIN_LIST;
    os_.write(buf, count);
    os_ << token::END_LIST;

    setState(os_.rdstate());

    return *this;
}


void Foam::OSstream::indent()
{
    for (unsigned short i = 0; i < indentLevel_*indentSize_; i++)
    {
        os_ << ' ';
    }
}


void Foam::OSstream::flush()
{
    os_.flush();
}


void Foam::OSstream::endl()
{
    write('\n');
    os_.flush();
}


std::ios_base::fmtflags Foam::OSstream::flags() const
{
    return os_.flags();
}


std::ios_base::fmtflags Foam::OSstream::flags(const ios_base::fmtflags f)
{
    return os_.flags(f);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int Foam::OSstream::width() const
{
    return os_.width();
}


int Foam::OSstream::width(const int w)
{
    return os_.width(w);
}


int Foam::OSstream::precision() const
{
    return os_.precision();
}


int Foam::OSstream::precision(const int p)
{
    return os_.precision(p);
}


// ************************************************************************* //
