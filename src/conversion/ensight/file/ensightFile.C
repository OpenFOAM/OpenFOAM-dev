/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "ensightFile.H"
#include <sstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::ensightFile::allowUndef_ = false;

Foam::scalar Foam::ensightFile::undefValue_ = Foam::floatScalarVGreat;

// default is width 8
Foam::string Foam::ensightFile::mask_ = "********";

Foam::string Foam::ensightFile::dirFmt_ = "%08d";


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::string Foam::ensightFile::mask()
{
    return mask_;
}


Foam::string Foam::ensightFile::subDir(const label n)
{
    char buf[32];

    sprintf(buf, dirFmt_.c_str(), n);
    return buf;
}


void Foam::ensightFile::subDirWidth(const label n)
{
    // enforce max limit to avoid buffer overflow in subDir()
    if (n < 1 || n > 31)
    {
        return;
    }

    // appropriate printf format
    std::ostringstream oss;
    oss << "%0" << n << "d";
    dirFmt_ = oss.str();

    // set mask accordingly
    mask_.resize(n, '*');
}


Foam::label Foam::ensightFile::subDirWidth()
{
    return mask_.size();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightFile::ensightFile
(
    const fileName& filePath,
    IOstream::streamFormat format
)
:
    OFstream(filePath, format)
{
    // ascii formatting specs
    setf
    (
        ios_base::scientific,
        ios_base::floatfield
    );
    precision(5);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightFile::~ensightFile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ensightFile::allowUndef()
{
    return allowUndef_;
}


bool Foam::ensightFile::allowUndef(bool value)
{
    bool old = allowUndef_;
    allowUndef_ = value;
    return old;
}


Foam::scalar Foam::ensightFile::undefValue(const scalar value)
{
    // enable its use too
    allowUndef_ = true;

    scalar old = undefValue_;
    undefValue_ = value;
    return old;
}


Foam::Ostream& Foam::ensightFile::write
(
    const char* buf,
    std::streamsize count
)
{
    stdStream().write(buf, count);
    return *this;
}


Foam::Ostream& Foam::ensightFile::write(const char* value)
{
    return write(string(value));
}


Foam::Ostream& Foam::ensightFile::write(const string& value)
{
    char buf[80];

    for (string::size_type i = 0; i < 80; ++i)
    {
        buf[i] = 0;
    }

    string::size_type n = value.size();
    if (n >= 80)
    {
        n = 79;
    }

    for (string::size_type i = 0; i < n; ++i)
    {
        buf[i] = value[i];
    }

    if (format() == IOstream::BINARY)
    {
        write
        (
            reinterpret_cast<char const *>(buf),
            sizeof(buf)
        );
    }
    else
    {
        stdStream() << buf;
    }

    return *this;
}


Foam::Ostream& Foam::ensightFile::write(const label value)
{
    if (format() == IOstream::BINARY)
    {
        unsigned int ivalue(value);

        write
        (
            reinterpret_cast<char const *>(&ivalue),
            sizeof(ivalue)
        );
    }
    else
    {
        stdStream().width(10);
        stdStream() << value;
    }

    return *this;
}


Foam::Ostream& Foam::ensightFile::write
(
    const label value,
    const label fieldWidth
)
{
    if (format() == IOstream::BINARY)
    {
        unsigned int ivalue(value);

        write
        (
            reinterpret_cast<char const *>(&ivalue),
            sizeof(ivalue)
        );
    }
    else
    {
        stdStream().width(fieldWidth);
        stdStream() << value;
    }

    return *this;
}


Foam::Ostream& Foam::ensightFile::write(const scalar value)
{
    float fvalue(value);

    if (format() == IOstream::BINARY)
    {
        write
        (
            reinterpret_cast<char const *>(&fvalue),
            sizeof(fvalue)
        );
    }
    else
    {
        stdStream().width(12);
        stdStream() << fvalue;
    }

    return *this;
}


void Foam::ensightFile::newline()
{
    if (format() == IOstream::ASCII)
    {
        stdStream() << nl;
    }
}


Foam::Ostream& Foam::ensightFile::writeUndef()
{
    write(undefValue_);
    return *this;
}


Foam::Ostream& Foam::ensightFile::writeKeyword(const string& key)
{
    if (allowUndef_)
    {
        write(string(key + " undef"));
        newline();
        write(undefValue_);
        newline();
    }
    else
    {
        write(key);
        newline();
    }
    return *this;
}


Foam::Ostream& Foam::ensightFile::writeBinaryHeader()
{
    if (format() == IOstream::BINARY)
    {
        write("C Binary");
    }

    return *this;
}


// ************************************************************************* //
