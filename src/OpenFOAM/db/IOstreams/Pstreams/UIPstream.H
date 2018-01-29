/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::UIPstream

Description
    Input inter-processor communications stream operating on external
    buffer.

SourceFiles
    UIPstream.C

\*---------------------------------------------------------------------------*/

#include "Pstream.H"

#ifndef UIPstream_H
#define UIPstream_H

#include "UPstream.H"
#include "Istream.H"
#include "PstreamBuffers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class UIPstream Declaration
\*---------------------------------------------------------------------------*/

class UIPstream
:
    public UPstream,
    public Istream
{
    // Private data

        int fromProcNo_;

        DynamicList<char>& externalBuf_;

        label& externalBufPosition_;

        const int tag_;

        const label comm_;

        const bool clearAtEnd_;

        int messageSize_;


    // Private Member Functions

        //- Check the bufferPosition against messageSize_ for EOF
        inline void checkEof();

        //- Read a T from the transfer buffer
        template<class T>
        inline void readFromBuffer(T&);

        //- Read data from the transfer buffer
        inline void readFromBuffer(void* data, size_t count, size_t align);


public:

    // Constructors

        //- Construct given process index to read from and optional buffer size,
        //  read format and IO version
        UIPstream
        (
            const commsTypes commsType,
            const int fromProcNo,
            DynamicList<char>& externalBuf,
            label& externalBufPosition,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            const bool clearAtEnd = false,   // destroy externalBuf if at end
            streamFormat format=BINARY,
            versionNumber version=currentVersion
        );

        //- Construct given buffers
        UIPstream(const int fromProcNo, PstreamBuffers&);


    //- Destructor
    ~UIPstream();


    // Member functions

        // Inquiry

            //- Return flags of output stream
            ios_base::fmtflags flags() const
            {
                return ios_base::fmtflags(0);
            }


        // Read functions

            //- Read into given buffer from given processor and return the
            //  message size
            static label read
            (
                const commsTypes commsType,
                const int fromProcNo,
                char* buf,
                const std::streamsize bufSize,
                const int tag = UPstream::msgType(),
                const label communicator = 0
            );

            //- Return next token from stream
            Istream& read(token&);

            //- Read a character
            Istream& read(char&);

            //- Read a word
            Istream& read(word&);

            // Read a string (including enclosing double-quotes)
            Istream& read(string&);

            //- Read a label
            Istream& read(label&);

            //- Read a floatScalar
            Istream& read(floatScalar&);

            //- Read a doubleScalar
            Istream& read(doubleScalar&);

            //- Read a longDoubleScalar
            Istream& read(longDoubleScalar&);

            //- Read binary block
            Istream& read(char*, std::streamsize);

            //- Rewind and return the stream so that it may be read again
            Istream& rewind();


        // Edit

            //- Set flags of stream
            ios_base::fmtflags flags(const ios_base::fmtflags)
            {
                return ios_base::fmtflags(0);
            }


        // Print

            //- Print description of IOstream to Ostream
            void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
