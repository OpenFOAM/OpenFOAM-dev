/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::ISstream

Description
    Generic input stream.

SourceFiles
    ISstreamI.H
    ISstream.C

\*---------------------------------------------------------------------------*/

#ifndef ISstream_H
#define ISstream_H

#include "Istream.H"
#include "fileName.H"
#include "DynamicList.H"
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class ISstream Declaration
\*---------------------------------------------------------------------------*/

class ISstream
:
    public Istream
{
    // Private Static Data

        //- Initial capacity of the character buffer
        static const int bufInitialCapacity = 1024;

        //- Amount of the buffer used in error messages. Less than the buffer
        //  length to facilitate readability.
        static const int bufErrorLength = 80;


    // Private Data

        //- The name of this stream
        fileName name_;

        //- Reference to the standard stream
        istream& is_;

        //- Character buffer
        DynamicList<char> buf_;


    // Private Member Functions

        char nextValid();

        //- Read a verbatim string (excluding block delimiters).
        Istream& readVerbatim(verbatimString&);

        //- Read a variable name into a string
        Istream& readVariable(string&);

        //- Read a delimited set of characters into a string
        Istream& readDelimited(string&, const char begin, const char end);

        //- Read a work token
        void readWordToken(token&);


public:

    // Constructors

        //- Construct as wrapper around istream
        inline ISstream
        (
            istream& is,
            const string& name,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED
        );


    //- Destructor
    virtual ~ISstream()
    {}


    // Member Functions

        // Inquiry

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


        // Read functions

            //- Low-level get character function.
            inline ISstream& get(char&);

            //- Low-level peek function.
            //  Does not remove the character from the stream.
            //  Returns the next character in the stream or EOF if the
            //  end of file is read.
            inline int peek();

            //- Low-level putback character function.
            inline ISstream& putback(const char&);

            //- Read line into a string
            //  with optional support for continuation lines
            ISstream& getLine(string&, const bool continuation = true);

            //- Read a '(...)' delimited set of characters into a string
            Istream& readList(string&);

            //- Read a '{...}' delimited set of characters into a string
            Istream& readBlock(string&);

            //- Return next token from stream
            virtual Istream& read(token&);

            //- Read a character
            virtual Istream& read(char&);

            //- Read a word
            virtual Istream& read(word&);

            //- Read a string (including enclosing double-quotes).
            //  Backslashes are retained, except when escaping double-quotes
            //  and an embedded newline character.
            virtual Istream& read(string&);

            //- Read a label
            virtual Istream& read(label&);

            //- Read a floatScalar
            virtual Istream& read(floatScalar&);

            //- Read a doubleScalar
            virtual Istream& read(doubleScalar&);

            //- Read a longDoubleScalar
            virtual Istream& read(longDoubleScalar&);

            //- Read binary block
            virtual Istream& read(char*, std::streamsize);

            //- Rewind and return the stream so that it may be read again
            virtual Istream& rewind();


        // Stream state functions

            //- Set flags of output stream
            virtual ios_base::fmtflags flags(const ios_base::fmtflags flags);


        // STL stream

            //- Access to underlying std::istream
            virtual istream& stdStream()
            {
                return is_;
            }

            //- Const access to underlying std::istream
            virtual const istream& stdStream() const
            {
                return is_;
            }


        // Print

            //- Print description of IOstream to Ostream
            virtual void print(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ISstream&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ISstreamI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
