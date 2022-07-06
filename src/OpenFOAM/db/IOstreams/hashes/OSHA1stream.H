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
    Foam::OSHA1stream

Description
    A Foam::OSstream for calculating SHA-1 digests

SourceFiles
    OSHA1stream.C

\*---------------------------------------------------------------------------*/

#ifndef OSHA1stream_H
#define OSHA1stream_H

#include "OSstream.H"
#include "SHA1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class sha1streambuf Declaration
\*---------------------------------------------------------------------------*/

class sha1streambuf
:
    public std::streambuf
{
    // Private Data

        //- The SHA-1 class which processes the strings
        SHA1 sha1_;


public:

    // Constructors

        //- Construct null
        sha1streambuf()
        {}


    // Member Functions

        //- Return the SHA-1
        SHA1& sha1()
        {
            return sha1_;
        }

        //- Writes characters from the array pointed to by s
        //  into the controlled output sequence, until either n characters
        //  have been written or the end of the output sequence is reached
        virtual std::streamsize xsputn(const char* s, std::streamsize n)
        {
            sha1_.append(s, n);
            return n;
        }

        //- Writes c to the current position of the put pointer (pptr),
        //  and advances that pointer one position forward.
        //  If c is EOF, no characters are written and the put pointer
        //  (pptr) is not advanced.
        //  Maybe used by xsputn (C++11).
        virtual int overflow(int c = EOF)
        {
            if (c != EOF)
            {
                const char cc = c;
                sha1_.append(&cc, 1);
            }

            return c;
        }
};


/*---------------------------------------------------------------------------*\
                         Class osha1stream Declaration
\*---------------------------------------------------------------------------*/

class osha1stream
:
    virtual public std::ios,
    public std::ostream
{
    // Private Data

        // SHA-1 stream buffer
        sha1streambuf sbuf_;


public:

    // Constructors

        //- Construct null
        osha1stream()
        :
            std::ostream(&sbuf_)
        {}


    // Member Functions

        //- Return the SHA-1
        SHA1& sha1()
        {
            return sbuf_.sha1();
        }

        //- Return the stream buffer
        sha1streambuf* rdbuf()
        {
            return &sbuf_;
        }
};


/*---------------------------------------------------------------------------*\
                         Class OSHA1stream Declaration
\*---------------------------------------------------------------------------*/

class OSHA1stream
:
    public OSstream
{

public:

    // Constructors

        //- Construct and set stream status
        OSHA1stream
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OSstream
            (
                *(new osha1stream),
                "OSHA1stream.sinkFile_",
                format,
                version
            )
        {}

        //- Disallow default bitwise copy construction
        OSHA1stream(const OSHA1stream&) = delete;


    //- Destructor
    ~OSHA1stream()
    {
        delete &dynamic_cast<osha1stream&>(stdStream());
    }


    // Member Functions

        //- Return the SHA-1
        SHA1& sha1()
        {
            return dynamic_cast<osha1stream&>(stdStream()).sha1();
        }

        //- Return the SHA-1 digest for the data processed until now
        SHA1Digest digest()
        {
            return sha1().digest();
        }

        //- Clear the SHA1 calculation
        void rewind()
        {
            sha1().clear();
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const OSHA1stream&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
