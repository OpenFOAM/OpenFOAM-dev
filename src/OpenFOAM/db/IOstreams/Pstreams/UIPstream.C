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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "UIPstream.H"
#include "int.H"
#include "token.H"

#include <cctype>


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::UIPstream::checkEof()
{
    if (externalBufPosition_ == messageSize_)
    {
        setEof();
    }
}


template<class T>
inline void Foam::UIPstream::readFromBuffer(T& t)
{
    const size_t align = sizeof(T);
    externalBufPosition_ = align + ((externalBufPosition_ - 1) & ~(align - 1));

    t = reinterpret_cast<T&>(externalBuf_[externalBufPosition_]);
    externalBufPosition_ += sizeof(T);
    checkEof();
}


inline void Foam::UIPstream::readFromBuffer
(
    void* data,
    size_t count,
    size_t align
)
{
    if (align > 1)
    {
        externalBufPosition_ =
            align
          + ((externalBufPosition_ - 1) & ~(align - 1));
    }

    const char* bufPtr = &externalBuf_[externalBufPosition_];
    char* dataPtr = reinterpret_cast<char*>(data);
    size_t i = count;
    while (i--) *dataPtr++ = *bufPtr++;
    externalBufPosition_ += count;
    checkEof();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::UIPstream::~UIPstream()
{
    if (clearAtEnd_ && eof())
    {
        if (debug)
        {
            Pout<< "UIPstream::~UIPstream() : tag:" << tag_
                << " fromProcNo:" << fromProcNo_
                << " clearing externalBuf_ of size "
                << externalBuf_.size()
                << " messageSize_:" << messageSize_ << endl;
        }
        externalBuf_.clearStorage();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Istream& Foam::UIPstream::read(token& t)
{
    // Return the put back token if it exists
    if (Istream::getBack(t))
    {
        return *this;
    }

    char c;

    // return on error
    if (!read(c))
    {
        t.setBad();
        return *this;
    }

    // Set the line number of this token to the current stream line number
    t.lineNumber() = lineNumber();

    // Analyse input starting with this character.
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
            t = token::punctuationToken(c);
            return *this;
        }

        // Word
        case token::WORD :
        {
            word* pval = new word;
            if (read(*pval))
            {
                if (token::compound::isCompound(*pval))
                {
                    t = token::compound::New(*pval, *this).ptr();
                    delete pval;
                }
                else
                {
                    t = pval;
                }
            }
            else
            {
                delete pval;
                t.setBad();
            }
            return *this;
        }

        // String
        case token::VERBATIMSTRING :
        {
            // Recurse to read actual string
            read(t);
            t.type() = token::VERBATIMSTRING;
            return *this;
        }
        case token::VARIABLE :
        {
            // Recurse to read actual string
            read(t);
            t.type() = token::VARIABLE;
            return *this;
        }
        case token::STRING :
        {
            string* pval = new string;
            if (read(*pval))
            {
                t = pval;
                if (c == token::VERBATIMSTRING)
                {
                    t.type() = token::VERBATIMSTRING;
                }
            }
            else
            {
                delete pval;
                t.setBad();
            }
            return *this;
        }

        // Label
        case token::LABEL :
        {
            label val;
            if (read(val))
            {
                t = val;
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // floatScalar
        case token::FLOAT_SCALAR :
        {
            floatScalar val;
            if (read(val))
            {
                t = val;
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // doubleScalar
        case token::DOUBLE_SCALAR :
        {
            doubleScalar val;
            if (read(val))
            {
                t = val;
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // longDoubleScalar
        case token::LONG_DOUBLE_SCALAR :
        {
            longDoubleScalar val;
            if (read(val))
            {
                t = val;
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // Character (returned as a single character word) or error
        default:
        {
            if (isalpha(c))
            {
                t = word(c);
                return *this;
            }

            setBad();
            t.setBad();

            return *this;
        }
    }
}


Foam::Istream& Foam::UIPstream::read(char& c)
{
    c = externalBuf_[externalBufPosition_];
    externalBufPosition_++;
    checkEof();
    return *this;
}


Foam::Istream& Foam::UIPstream::read(word& str)
{
    size_t len;
    readFromBuffer(len);
    str = &externalBuf_[externalBufPosition_];
    externalBufPosition_ += len + 1;
    checkEof();
    return *this;
}


Foam::Istream& Foam::UIPstream::read(string& str)
{
    size_t len;
    readFromBuffer(len);
    str = &externalBuf_[externalBufPosition_];
    externalBufPosition_ += len + 1;
    checkEof();
    return *this;
}


Foam::Istream& Foam::UIPstream::read(label& val)
{
    readFromBuffer(val);
    return *this;
}


Foam::Istream& Foam::UIPstream::read(floatScalar& val)
{
    readFromBuffer(val);
    return *this;
}


Foam::Istream& Foam::UIPstream::read(doubleScalar& val)
{
    readFromBuffer(val);
    return *this;
}


Foam::Istream& Foam::UIPstream::read(longDoubleScalar& val)
{
    readFromBuffer(val);
    return *this;
}


Foam::Istream& Foam::UIPstream::read(char* data, std::streamsize count)
{
    if (format() != BINARY)
    {
        FatalErrorInFunction
            << "stream format not binary"
            << Foam::abort(FatalError);
    }

    readFromBuffer(data, count, 8);
    return *this;
}


Foam::Istream& Foam::UIPstream::rewind()
{
    externalBufPosition_ = 0;
    return *this;
}


void Foam::UIPstream::print(Ostream& os) const
{
    os  << "Reading from processor " << fromProcNo_
        << " using communicator " << comm_
        <<  " and tag " << tag_
        << Foam::endl;
}


// ************************************************************************* //
