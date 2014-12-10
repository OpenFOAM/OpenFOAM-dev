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

\*---------------------------------------------------------------------------*/

#include "ISstream.H"
#include "int.H"
#include "token.H"
#include <cctype>

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

char Foam::ISstream::nextValid()
{
    char c = 0;

    while (true)
    {
        // Get next non-whitespace character
        while (get(c) && isspace(c))
        {}

        // Return if stream is bad - ie, previous get() failed
        if (bad() || isspace(c))
        {
            break;
        }

        // Is this the start of a C/C++ comment?
        if (c == '/')
        {
            if (!get(c))
            {
                // cannot get another character - return this one
                return '/';
            }

            if (c == '/')
            {
                // C++ style single-line comment - skip through past end-of-line
                while (get(c) && c != '\n')
                {}
            }
            else if (c == '*')
            {
                // within a C-style comment
                while (true)
                {
                    // search for end of C-style comment - '*/'
                    if (get(c) && c == '*')
                    {
                        if (get(c))
                        {
                            if (c == '/')
                            {
                                // matched '*/'
                                break;
                            }
                            else if (c == '*')
                            {
                                // check again
                                putback(c);
                            }
                        }
                    }

                    if (!good())
                    {
                        return 0;
                    }
                }
            }
            else
            {
                // The '/' did not start a C/C++ comment - return it
                putback(c);
                return '/';
            }
        }
        else
        {
            // a valid character - return it
            return c;
        }
    }

    return 0;
}


void Foam::ISstream::readWordToken(token& t)
{
    word* wPtr = new word;

    if (read(*wPtr).bad())
    {
        delete wPtr;
        t.setBad();
    }
    else if (token::compound::isCompound(*wPtr))
    {
        t = token::compound::New(*wPtr, *this).ptr();
        delete wPtr;
    }
    else
    {
        t = wPtr;
    }
}


Foam::Istream& Foam::ISstream::read(token& t)
{
    static const int maxLen = 128;
    static char buf[maxLen];

    // Return the put back token if it exists
    if (Istream::getBack(t))
    {
        return *this;
    }

    // Assume that the streams supplied are in working order.
    // Lines are counted by '\n'

    // Get next 'valid character': i.e. proceed through any whitespace
    // and/or comments until a semantically valid character is found

    char c = nextValid();

    // Set the line number of this token to the current stream line number
    t.lineNumber() = lineNumber();

    // return on error
    if (!c)
    {
        t.setBad();
        return *this;
    }

    // Analyse input starting with this character.
    switch (c)
    {
        // Check for punctuation first

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
        // NB: token::SUBTRACT handled later as the possible start of a Number
        case token::MULTIPLY :
        case token::DIVIDE :
        {
            t = token::punctuationToken(c);
            return *this;
        }


        // String: enclosed by double quotes.
        case token::BEGIN_STRING :
        {
            putback(c);
            string* sPtr = new string;

            if (read(*sPtr).bad())
            {
                delete sPtr;
                t.setBad();
            }
            else
            {
                t = sPtr;
            }

            return *this;
        }
        // Possible verbatim string or dictionary functionEntry
        case token::HASH :
        {
            char nextC;
            if (read(nextC).bad())
            {
                // Return hash as word
                t = token(word(c));
                return *this;
            }
            else if (nextC == token::BEGIN_BLOCK)
            {
                // Verbatim string
                string* sPtr = new string;

                if (readVerbatim(*sPtr).bad())
                {
                    delete sPtr;
                    t.setBad();
                }
                else
                {
                    t = sPtr;
                    t.type() = token::VERBATIMSTRING;
                }

                return *this;
            }
            else
            {
                // Word beginning with #
                putback(nextC);
                putback(c);

                readWordToken(t);

                return *this;
            }
        }

        case '$':
        {
            // Look ahead
            char nextC;
            if (read(nextC).bad())
            {
                // Return $ as word
                t = token(word(c));
                return *this;
            }
            else if (nextC == token::BEGIN_BLOCK)
            {
                putback(nextC);
                putback(c);

                string* sPtr = new string;

                if (readVariable(*sPtr).bad())
                {
                    delete sPtr;
                    t.setBad();
                }
                else
                {
                    t = sPtr;
                    t.type() = token::VARIABLE;
                }
                return *this;
            }
            else
            {
                putback(nextC);
                putback(c);
                readWordToken(t);
                return *this;
            }
        }

        // Number: integer or floating point
        //
        // ideally match the equivalent of this regular expression
        //
        //    /[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([Ee][-+]?[0-9]+)?/
        //
        case '-' :
        case '.' :
        case '0' : case '1' : case '2' : case '3' : case '4' :
        case '5' : case '6' : case '7' : case '8' : case '9' :
        {
            bool asLabel = (c != '.');

            int nChar = 0;
            buf[nChar++] = c;

            // get everything that could resemble a number and let
            // readScalar determine the validity
            while
            (
                is_.get(c)
             && (
                    isdigit(c)
                 || c == '+'
                 || c == '-'
                 || c == '.'
                 || c == 'E'
                 || c == 'e'
                )
            )
            {
                if (asLabel)
                {
                    asLabel = isdigit(c);
                }

                buf[nChar++] = c;
                if (nChar == maxLen)
                {
                    // runaway argument - avoid buffer overflow
                    buf[maxLen-1] = '\0';

                    FatalIOErrorIn("ISstream::read(token&)", *this)
                        << "number '" << buf << "...'\n"
                        << "    is too long (max. " << maxLen << " characters)"
                        << exit(FatalIOError);

                    t.setBad();
                    return *this;
                }
            }
            buf[nChar] = '\0';

            setState(is_.rdstate());
            if (is_.bad())
            {
                t.setBad();
            }
            else
            {
                is_.putback(c);

                if (nChar == 1 && buf[0] == '-')
                {
                    // a single '-' is punctuation
                    t = token::punctuationToken(token::SUBTRACT);
                }
                else
                {
                    if (asLabel)
                    {
                        label labelVal;
                        if (readLabel(buf, labelVal))
                        {
                            t = labelVal;
                        }
                        else
                        {
                            // Maybe too big? Try as scalar
                            scalar scalarVal;
                            if (readScalar(buf, scalarVal))
                            {
                                t = scalarVal;
                            }
                            else
                            {
                                t.setBad();
                            }
// ---------------------------------------
// this would also be possible if desired:
// ---------------------------------------
//                        // return as a label when possible
//                        // eg, 1E6 -> 1000000
//                        if (scalarVal <= labelMax && scalarVal >= labelMin)
//                        {
//                            label labelVal(scalarVal);
//
//                            if (labelVal == scalarVal)
//                            {
//                                t = labelVal;
//                            }
//                        }

                        }
                    }
                    else
                    {
                        scalar scalarVal;
                        if (readScalar(buf, scalarVal))
                        {
                            t = scalarVal;
                        }
                        else
                        {
                            t.setBad();
                        }
                    }
                }
            }

            return *this;
        }


        // Should be a word (which can also be a single character)
        default:
        {
            putback(c);
            readWordToken(t);

            return *this;
        }
    }
}


Foam::Istream& Foam::ISstream::read(char& c)
{
    c = nextValid();
    return *this;
}


Foam::Istream& Foam::ISstream::read(word& str)
{
    static const int maxLen = 1024;
    static const int errLen = 80; // truncate error message for readability
    static char buf[maxLen];

    register int nChar = 0;
    register int listDepth = 0;
    char c;

    while (get(c) && word::valid(c))
    {
        if (c == token::BEGIN_LIST)
        {
            listDepth++;
        }
        else if (c == token::END_LIST)
        {
            if (listDepth)
            {
                listDepth--;
            }
            else
            {
                break;
            }
        }

        buf[nChar++] = c;
        if (nChar == maxLen)
        {
            buf[errLen] = '\0';

            FatalIOErrorIn("ISstream::read(word&)", *this)
                << "word '" << buf << "...'\n"
                << "    is too long (max. " << maxLen << " characters)"
                << exit(FatalIOError);

            return *this;
        }
    }

    // we could probably skip this check
    if (bad())
    {
        buf[errLen] = buf[nChar] = '\0';

        FatalIOErrorIn("ISstream::read(word&)", *this)
            << "problem while reading word '" << buf << "...' after "
            << nChar << " characters\n"
            << exit(FatalIOError);

        return *this;
    }

    if (nChar == 0)
    {
        FatalIOErrorIn("ISstream::read(word&)", *this)
            << "invalid first character found : " << c
            << exit(FatalIOError);
    }

    // done reading
    buf[nChar] = '\0';
    str = buf;
    putback(c);

    return *this;
}


Foam::Istream& Foam::ISstream::read(string& str)
{
    static const int maxLen = 1024;
    static const int errLen = 80; // truncate error message for readability
    static char buf[maxLen];

    char c;

    if (!get(c))
    {
        FatalIOErrorIn("ISstream::read(string&)", *this)
            << "cannot read start of string"
            << exit(FatalIOError);

        return *this;
    }

    // Note, we could also handle single-quoted strings here (if desired)
    if (c != token::BEGIN_STRING)
    {
        FatalIOErrorIn("ISstream::read(string&)", *this)
            << "Incorrect start of string character found : " << c
            << exit(FatalIOError);

        return *this;
    }

    register int nChar = 0;
    bool escaped = false;

    while (get(c))
    {
        if (c == token::END_STRING)
        {
            if (escaped)
            {
                escaped = false;
                nChar--;    // overwrite backslash
            }
            else
            {
                // done reading
                buf[nChar] = '\0';
                str = buf;
                return *this;
            }
        }
        else if (c == token::NL)
        {
            if (escaped)
            {
                escaped = false;
                nChar--;    // overwrite backslash
            }
            else
            {
                buf[errLen] = buf[nChar] = '\0';

                FatalIOErrorIn("ISstream::read(string&)", *this)
                    << "found '\\n' while reading string \""
                    << buf << "...\""
                    << exit(FatalIOError);

                return *this;
            }
        }
        else if (c == '\\')
        {
            escaped = !escaped;    // toggle state (retains backslashes)
        }
        else
        {
            escaped = false;
        }

        buf[nChar++] = c;
        if (nChar == maxLen)
        {
            buf[errLen] = '\0';

            FatalIOErrorIn("ISstream::read(string&)", *this)
                << "string \"" << buf << "...\"\n"
                << "    is too long (max. " << maxLen << " characters)"
                << exit(FatalIOError);

            return *this;
        }
    }


    // don't worry about a dangling backslash if string terminated prematurely
    buf[errLen] = buf[nChar] = '\0';

    FatalIOErrorIn("ISstream::read(string&)", *this)
        << "problem while reading string \"" << buf << "...\""
        << exit(FatalIOError);

    return *this;
}


// Special handling of '{' in variables
Foam::Istream& Foam::ISstream::readVariable(string& str)
{
    static const int maxLen = 1024;
    static const int errLen = 80; // truncate error message for readability
    static char buf[maxLen];

    register int nChar = 0;
    register int blockCount = 0;
    char c;

    if (!get(c) || c != '$')
    {
        FatalIOErrorIn("ISstream::readVariable(string&)", *this)
            << "invalid first character found : " << c
            << exit(FatalIOError);
    }

    buf[nChar++] = c;

    // Read next character to see if '{'
    if (get(c) && c == token::BEGIN_BLOCK)
    {
        // Read, counting brackets
        buf[nChar++] = c;

        while
        (
            get(c)
         && (
                c == token::BEGIN_BLOCK
             || c == token::END_BLOCK
             || word::valid(c)
            )
        )
        {
            buf[nChar++] = c;
            if (nChar == maxLen)
            {
                buf[errLen] = '\0';

                FatalIOErrorIn("ISstream::readVariable(string&)", *this)
                    << "word '" << buf << "...'\n"
                    << "    is too long (max. " << maxLen << " characters)"
                    << exit(FatalIOError);

                return *this;
            }

            if (c == token::BEGIN_BLOCK)
            {
                blockCount++;
            }
            else if (c == token::END_BLOCK)
            {
                if (blockCount)
                {
                    blockCount--;
                }
                else
                {
                    break;
                }
            }
        }
    }
    else
    {
        buf[nChar++] = c;

        while (get(c) && word::valid(c))
        {
            buf[nChar++] = c;
            if (nChar == maxLen)
            {
                buf[errLen] = '\0';

                FatalIOErrorIn("ISstream::readVariable(string&)", *this)
                    << "word '" << buf << "...'\n"
                    << "    is too long (max. " << maxLen << " characters)"
                    << exit(FatalIOError);

                return *this;
            }
        }
    }

    // we could probably skip this check
    if (bad())
    {
        buf[errLen] = buf[nChar] = '\0';

        FatalIOErrorIn("ISstream::readVariable(string&)", *this)
            << "problem while reading string '" << buf << "...' after "
            << nChar << " characters\n"
            << exit(FatalIOError);

        return *this;
    }

    if (nChar == 0)
    {
        FatalIOErrorIn("ISstream::readVariable(string&)", *this)
            << "invalid first character found : " << c
            << exit(FatalIOError);
    }

    // done reading
    buf[nChar] = '\0';
    str = buf;

    // Note: check if we exited due to '}' or just !word::valid.
    if (c != token::END_BLOCK)
    {
        putback(c);
    }

    return *this;
}


Foam::Istream& Foam::ISstream::readVerbatim(string& str)
{
    static const int maxLen = 8000;
    static const int errLen = 80; // truncate error message for readability
    static char buf[maxLen];

    char c;

    register int nChar = 0;

    while (get(c))
    {
        if (c == token::HASH)
        {
            char nextC;
            get(nextC);
            if (nextC == token::END_BLOCK)
            {
                buf[nChar] = '\0';
                str = buf;
                return *this;
            }
            else
            {
                putback(nextC);
            }
        }

        buf[nChar++] = c;
        if (nChar == maxLen)
        {
            buf[errLen] = '\0';

            FatalIOErrorIn("ISstream::readVerbatim(string&)", *this)
                << "string \"" << buf << "...\"\n"
                << "    is too long (max. " << maxLen << " characters)"
                << exit(FatalIOError);

            return *this;
        }
    }


    // don't worry about a dangling backslash if string terminated prematurely
    buf[errLen] = buf[nChar] = '\0';

    FatalIOErrorIn("ISstream::readVerbatim(string&)", *this)
        << "problem while reading string \"" << buf << "...\""
        << exit(FatalIOError);

    return *this;
}


Foam::Istream& Foam::ISstream::read(label& val)
{
    is_ >> val;
    setState(is_.rdstate());
    return *this;
}


Foam::Istream& Foam::ISstream::read(floatScalar& val)
{
    is_ >> val;
    setState(is_.rdstate());
    return *this;
}


Foam::Istream& Foam::ISstream::read(doubleScalar& val)
{
    is_ >> val;
    setState(is_.rdstate());
    return *this;
}


// read binary block
Foam::Istream& Foam::ISstream::read(char* buf, std::streamsize count)
{
    if (format() != BINARY)
    {
        FatalIOErrorIn("ISstream::read(char*, std::streamsize)", *this)
            << "stream format not binary"
            << exit(FatalIOError);
    }

    readBegin("binaryBlock");
    is_.read(buf, count);
    readEnd("binaryBlock");

    setState(is_.rdstate());

    return *this;
}


Foam::Istream& Foam::ISstream::rewind()
{
    stdStream().rdbuf()->pubseekpos(0);

    return *this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


std::ios_base::fmtflags Foam::ISstream::flags() const
{
    return is_.flags();
}


std::ios_base::fmtflags Foam::ISstream::flags(const ios_base::fmtflags f)
{
    return is_.flags(f);
}


// ************************************************************************* //
