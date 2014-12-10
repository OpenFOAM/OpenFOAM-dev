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

#include "IOstream.H"
#include "error.H"
#include "Switch.H"
#include <sstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::fileName Foam::IOstream::name_("IOstream");


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

Foam::IOstream::streamFormat
Foam::IOstream::formatEnum(const word& format)
{
    if (format == "ascii")
    {
        return IOstream::ASCII;
    }
    else if (format == "binary")
    {
        return IOstream::BINARY;
    }
    else
    {
        WarningIn("IOstream::formatEnum(const word&)")
            << "bad format specifier '" << format << "', using 'ascii'"
            << endl;

        return IOstream::ASCII;
    }
}


Foam::IOstream::compressionType
Foam::IOstream::compressionEnum(const word& compression)
{
    // get Switch (bool) value, but allow it to fail
    Switch sw(compression, true);

    if (sw.valid())
    {
        return sw ? IOstream::COMPRESSED : IOstream::UNCOMPRESSED;
    }
    else if (compression == "uncompressed")
    {
        return IOstream::UNCOMPRESSED;
    }
    else if (compression == "compressed")
    {
        return IOstream::COMPRESSED;
    }
    else
    {
        WarningIn("IOstream::compressionEnum(const word&)")
            << "bad compression specifier '" << compression
            << "', using 'uncompressed'"
            << endl;

        return IOstream::UNCOMPRESSED;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::IOstream::check(const char* operation) const
{
    if (bad())
    {
        FatalIOErrorIn
        (
            "IOstream::check(const char*) const", *this
        )   << "error in IOstream " << name() << " for operation " << operation
            << exit(FatalIOError);
    }

    return !bad();
}


void Foam::IOstream::fatalCheck(const char* operation) const
{
    if (bad())
    {
        FatalIOErrorIn
        (
            "IOstream::fatalCheck(const char*) const", *this
        )   << "error in IOstream " << name() << " for operation " << operation
            << exit(FatalIOError);
    }
}


Foam::string Foam::IOstream::versionNumber::str() const
{
    std::ostringstream os;
    os.precision(1);
    os.setf(ios_base::fixed, ios_base::floatfield);
    os  << versionNumber_;
    return os.str();
}


void Foam::IOstream::print(Ostream& os) const
{
    os  << "IOstream: " << "Version "  << version_ << ", format ";

    switch (format_)
    {
        case ASCII:
            os  << "ASCII";
        break;

        case BINARY:
            os  << "BINARY";
        break;
    }

    os  << ", line "       << lineNumber();

    if (opened())
    {
        os  << ", OPENED";
    }

    if (closed())
    {
        os  << ", CLOSED";
    }

    if (good())
    {
        os  << ", GOOD";
    }

    if (eof())
    {
        os  << ", EOF";
    }

    if (fail())
    {
        os  << ", FAIL";
    }

    if (bad())
    {
        os  << ", BAD";
    }

    os  << endl;
}


void Foam::IOstream::print(Ostream& os, const int streamState) const
{
    if (streamState == ios_base::goodbit)
    {
        os  << "ios_base::goodbit set : the last operation on stream succeeded"
            << endl;
    }
    else if (streamState & ios_base::badbit)
    {
        os  << "ios_base::badbit set : characters possibly lost"
            << endl;
    }
    else if (streamState & ios_base::failbit)
    {
        os  << "ios_base::failbit set : some type of formatting error"
            << endl;
    }
    else if (streamState & ios_base::eofbit)
    {
        os  << "ios_base::eofbit set : at end of stream"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const IOstream::streamFormat& sf)
{
    if (sf == IOstream::ASCII)
    {
        os  << "ascii";
    }
    else
    {
        os  << "binary";
    }

    return os;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const IOstream::versionNumber& vn)
{
    os  << vn.str().c_str();
    return os;
}



namespace Foam
{

#  if defined (__GNUC__)
   template<>
#  endif
   Ostream& operator<<(Ostream& os, const InfoProxy<IOstream>& ip)
   {
       ip.t_.print(os);
       return os;
   }

}  // end namespace


// ************************************************************************* //
