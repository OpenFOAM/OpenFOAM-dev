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
    Write primitive and binary block from UOPstream

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "UOPstream.H"
#include "int.H"
#include "token.H"

#include <cctype>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
inline void Foam::UOPstream::writeToBuffer(const T& t)
{
    writeToBuffer(&t, sizeof(T), sizeof(T));
}


inline void Foam::UOPstream::writeToBuffer(const char& c)
{
    if (!sendBuf_.capacity())
    {
        sendBuf_.setCapacity(1000);
    }
    sendBuf_.append(c);
}


inline void Foam::UOPstream::writeToBuffer
(
    const void* data,
    size_t count,
    size_t align
)
{
    if (!sendBuf_.capacity())
    {
        sendBuf_.setCapacity(1000);
    }

    label alignedPos = sendBuf_.size();

    if (align > 1)
    {
        // Align bufPosition. Pads sendBuf_.size() - oldPos characters.
        alignedPos = align + ((sendBuf_.size() - 1) & ~(align - 1));
    }

    // Extend if necessary
    sendBuf_.setSize(alignedPos + count);

    register const char* dataPtr = reinterpret_cast<const char*>(data);
    register size_t i = count;
    while (i--) sendBuf_[alignedPos++] = *dataPtr++;
}



// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::UOPstream::UOPstream
(
    const commsTypes commsType,
    const int toProcNo,
    DynamicList<char>& sendBuf,
    const int tag,
    const label comm,
    const bool sendAtDestruct,
    streamFormat format,
    versionNumber version
)
:
    UPstream(commsType),
    Ostream(format, version),
    toProcNo_(toProcNo),
    sendBuf_(sendBuf),
    tag_(tag),
    comm_(comm),
    sendAtDestruct_(sendAtDestruct)
{
    setOpened();
    setGood();
}


Foam::UOPstream::UOPstream(const int toProcNo, PstreamBuffers& buffers)
:
    UPstream(buffers.commsType_),
    Ostream(buffers.format_, buffers.version_),
    toProcNo_(toProcNo),
    sendBuf_(buffers.sendBuf_[toProcNo]),
    tag_(buffers.tag_),
    comm_(buffers.comm_),
    sendAtDestruct_(buffers.commsType_ != UPstream::nonBlocking)
{
    setOpened();
    setGood();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::UOPstream::~UOPstream()
{
    if (sendAtDestruct_)
    {
        if
        (
           !UOPstream::write
            (
                commsType_,
                toProcNo_,
                sendBuf_.begin(),
                sendBuf_.size(),
                tag_,
                comm_
            )
        )
        {
            FatalErrorIn("UOPstream::~UOPstream()")
                << "Failed sending outgoing message of size " << sendBuf_.size()
                << " to processor " << toProcNo_
                << Foam::abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::UOPstream::write(const token& t)
{
    // Raw token output only supported for verbatim strings for now
    if (t.type() == token::VERBATIMSTRING)
    {
        write(char(token::VERBATIMSTRING));
        write(t.stringToken());
    }
    else if (t.type() == token::VARIABLE)
    {
        write(char(token::VARIABLE));
        write(t.stringToken());
    }
    else
    {
        notImplemented("Ostream& UOPstream::write(const token&)");
        setBad();
    }
    return *this;
}


Foam::Ostream& Foam::UOPstream::write(const char c)
{
    if (!isspace(c))
    {
        writeToBuffer(c);
    }

    return *this;
}


Foam::Ostream& Foam::UOPstream::write(const char* str)
{
    word nonWhiteChars(string::validate<word>(str));

    if (nonWhiteChars.size() == 1)
    {
        return write(nonWhiteChars.c_str()[1]);
    }
    else if (nonWhiteChars.size())
    {
        return write(nonWhiteChars);
    }
    else
    {
        return *this;
    }
}


Foam::Ostream& Foam::UOPstream::write(const word& str)
{
    write(char(token::WORD));

    size_t len = str.size();
    writeToBuffer(len);
    writeToBuffer(str.c_str(), len + 1, 1);

    return *this;
}


Foam::Ostream& Foam::UOPstream::write(const string& str)
{
    write(char(token::STRING));

    size_t len = str.size();
    writeToBuffer(len);
    writeToBuffer(str.c_str(), len + 1, 1);

    return *this;
}


Foam::Ostream& Foam::UOPstream::writeQuoted
(
    const std::string& str,
    const bool quoted
)
{
    if (quoted)
    {
        write(char(token::STRING));
    }
    else
    {
        write(char(token::WORD));
    }

    size_t len = str.size();
    writeToBuffer(len);
    writeToBuffer(str.c_str(), len + 1, 1);

    return *this;
}


Foam::Ostream& Foam::UOPstream::write(const label val)
{
    write(char(token::LABEL));
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstream::write(const floatScalar val)
{
    write(char(token::FLOAT_SCALAR));
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstream::write(const doubleScalar val)
{
    write(char(token::DOUBLE_SCALAR));
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstream::write(const char* data, std::streamsize count)
{
    if (format() != BINARY)
    {
        FatalErrorIn("Ostream::write(const char*, std::streamsize)")
            << "stream format not binary"
            << Foam::abort(FatalError);
    }

    writeToBuffer(data, count, 8);

    return *this;
}


void Foam::UOPstream::print(Ostream& os) const
{
    os  << "Writing from processor " << toProcNo_
        << " to processor " << myProcNo() << " in communicator " << comm_
        << " and tag " << tag_ << Foam::endl;
}


// ************************************************************************* //
