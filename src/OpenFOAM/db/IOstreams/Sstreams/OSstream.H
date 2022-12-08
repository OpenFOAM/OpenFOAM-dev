/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

Class
    Foam::OSstream

Description
    Generic output stream.

SourceFiles
    OSstreamI.H
    OSstream.C

\*---------------------------------------------------------------------------*/

#ifndef OSstream_H
#define OSstream_H

#include "Ostream.H"
#include "fileName.H"
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class OSstream Declaration
\*---------------------------------------------------------------------------*/

class OSstream
:
    public Ostream
{
    // Private Data

        fileName name_;
        ostream& os_;


public:

    // Constructors

        //- Set stream status
        OSstream
        (
            ostream& os,
            const string& name,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED
        );

        //- Copy construct
        OSstream(const OSstream&) = default;


    // Member Functions

        // Enquiry

            //- Return the name of the stream
            //  Useful for Fstream to return the filename
            virtual const fileName& name() const
            {
                return name_;
            }

            //- Return non-const access to the name of the stream
            //  Useful to alter the stream name
            virtual fileName& name()
            {
                return name_;
            }

            //- Return flags of output stream
            virtual ios_base::fmtflags flags() const;


        // Write functions

            //- Write character
            virtual Ostream& write(const char);

            //- Write character string
            virtual Ostream& write(const char*);

            //- Write word
            virtual Ostream& write(const word&);

            //- Write string with double quotes
            virtual Ostream& write(const string&);

            //- Write verbatimString with #{ }#
            virtual Ostream& write(const verbatimString&);

            //- Write std::string with optional double quotes.
            virtual Ostream& writeQuoted
            (
                const std::string&,
                const bool quoted=true
            );

            //- Write int32_t
            virtual Ostream& write(const int32_t);

            //- Write int64_t
            virtual Ostream& write(const int64_t);

            //- Write floatScalar
            virtual Ostream& write(const floatScalar);

            //- Write doubleScalar
            virtual Ostream& write(const doubleScalar);

            //- Write longDoubleScalar
            virtual Ostream& write(const longDoubleScalar);

            //- Write binary block
            virtual Ostream& write(const char*, std::streamsize);

            //- Add indentation characters
            virtual void indent();


        // Stream state functions

            //- Set flags of output stream
            virtual ios_base::fmtflags flags(const ios_base::fmtflags flags);

            //- Flush stream
            virtual void flush();

            //- Add newline and flush stream
            virtual void endl();

            //- Get width of output field
            virtual int width() const;

            //- Set width of output field (and return old width)
            virtual int width(const int);

            //- Get precision of output field
            virtual int precision() const;

            //- Set precision of output field (and return old precision)
            virtual int precision(const int);


        // STL stream

            //- Access to underlying std::ostream
            virtual ostream& stdStream()
            {
                return os_;
            }

            //- Const access to underlying std::ostream
            virtual const ostream& stdStream() const
            {
                return os_;
            }


        // Print

            //- Print description of IOstream to Ostream
            virtual void print(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const OSstream&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "OSstreamI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
