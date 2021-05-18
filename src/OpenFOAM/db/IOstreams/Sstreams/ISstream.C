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

\*---------------------------------------------------------------------------*/

#include "ISstream.H"
#include "int.H"
#include "token.H"
#include "DynamicList.H"
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
                // Cannot get another character - return this one
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
                // Within a C-style comment
                while (true)
                {
                    // Search for end of C-style comment - '*/'
                    if (get(c) && c == '*')
                    {
                        if (get(c))
                        {
                            if (c == '/')
                            {
                                // Matched '*/'
                                break;
                            }
                            else if (c == '*')
                            {
                                // Check again
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
            // A valid character - return it
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

    // Return on error
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
                verbatimString* vsPtr = new verbatimString;

                if (readVerbatim(*vsPtr).bad())
                {
                    delete vsPtr;
                    t.setBad();
                }
                else
                {
                    t = vsPtr;
                }
                return *this;
            }
            else
            {
                // Function name beginning with a #
                putback(nextC);
                putback(c);

                functionName* fnPtr = new functionName;

                if (read(*fnPtr).bad())
                {
                    delete fnPtr;
                    t.setBad();
                }
                else
                {
                    t = fnPtr;
                }
                return *this;
            }
        }

        case '$' :
        {
            // Look ahead
            char nextC;
            if (read(nextC).bad())
            {
                // Return $ as a word
                t = token(word(c));
                return *this;
            }
            else
            {
                putback(nextC);
                putback(c);

                variable* vPtr = new variable;

                if (readVariable(*vPtr).bad())
                {
                    delete vPtr;
                    t.setBad();
                }
                else
                {
                    t = vPtr;
                }

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

            buf_.clear();
            buf_.append(c);

            // Get everything that could resemble a number and let
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

                buf_.append(c);
            }

            buf_.append('\0');

            setState(is_.rdstate());
            if (is_.bad())
            {
                t.setBad();
            }
            else
            {
                is_.putback(c);

                if (buf_.size() == 2 && buf_[0] == '-')
                {
                    // A single '-' is punctuation
                    t = token::punctuationToken(token::SUBTRACT);
                }
                else
                {
                    if (asLabel)
                    {
                        label labelVal = 0;
                        if (Foam::read(buf_.cdata(), labelVal))
                        {
                            t = labelVal;
                        }
                        else
                        {
                            // Maybe too big? Try as scalar
                            scalar scalarVal;
                            if (readScalar(buf_.cdata(), scalarVal))
                            {
                                t = scalarVal;
                            }
                            else
                            {
                                t.setBad();
                            }
                        }
                    }
                    else
                    {
                        scalar scalarVal;
                        if (readScalar(buf_.cdata(), scalarVal))
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
    buf_.clear();

    int listDepth = 0;
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

        buf_.append(c);
    }

    if (bad())
    {
        // Truncate the string to the maximum error length
        buf_.data()[bufErrorLength] = buf_.last() = '\0';

        FatalIOErrorInFunction(*this)
            << "problem while reading word '" << buf_.cdata() << "...' after "
            << buf_.size() << " characters\n"
            << exit(FatalIOError);

        return *this;
    }

    if (buf_.empty())
    {
        FatalIOErrorInFunction(*this)
            << "invalid first character found : " << c
            << exit(FatalIOError);
    }

    // Append end string character
    buf_.append('\0');
    str = buf_.cdata();
    putback(c);

    return *this;
}


Foam::Istream& Foam::ISstream::read(string& str)
{
    buf_.clear();

    char c;

    if (!get(c))
    {
        FatalIOErrorInFunction(*this)
            << "cannot read start of string"
            << exit(FatalIOError);

        return *this;
    }

    // Note, we could also handle single-quoted strings here (if desired)
    if (c != token::BEGIN_STRING)
    {
        FatalIOErrorInFunction(*this)
            << "Incorrect start of string character found : " << c
            << exit(FatalIOError);

        return *this;
    }

    bool escaped = false;

    while (get(c))
    {
        if (c == token::END_STRING)
        {
            if (escaped)
            {
                escaped = false;
                buf_.remove(); // Overwrite backslash
            }
            else
            {
                // Append end string character
                buf_.append('\0');
                str = buf_.cdata();
                return *this;
            }
        }
        else if (c == token::NL)
        {
            if (escaped)
            {
                escaped = false;
                buf_.remove(); // Overwrite backslash
            }
            else
            {
                // Truncate the string to the maximum error length
                buf_.data()[bufErrorLength] = buf_.last() = '\0';

                FatalIOErrorInFunction(*this)
                    << "found '\\n' while reading string \""
                    << buf_.cdata() << "...\""
                    << exit(FatalIOError);

                return *this;
            }
        }
        else if (c == '\\')
        {
            escaped = !escaped;    // Toggle state (retains backslashes)
        }
        else
        {
            escaped = false;
        }

        buf_.append(c);
    }

    // Truncate the string to the maximum error length
    buf_.data()[bufErrorLength] = buf_.last() = '\0';

    FatalIOErrorInFunction(*this)
        << "problem while reading string \"" << buf_.cdata() << "...\""
        << exit(FatalIOError);

    return *this;
}


Foam::Istream& Foam::ISstream::readVariable(string& str)
{
    buf_.clear();

    char c;

    if (!get(c) || c != '$')
    {
        FatalIOErrorInFunction(*this)
            << "invalid first character found : " << c
            << exit(FatalIOError);
    }

    buf_.append(c);

    // Read next character to see if '{'
    if (get(c) && c == token::BEGIN_BLOCK)
    {
        // Read, counting brackets
        buf_.append(c);

        while
        (
            get(c)
         && (
                c == token::BEGIN_BLOCK
             || c == token::END_BLOCK
             || variable::valid(c)
            )
        )
        {
            buf_.append(c);

            int blockCount = 0;

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
        buf_.append(c);

        int listDepth = 0;

        while (get(c) && variable::valid(c))
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

            buf_.append(c);
        }
    }

    if (bad())
    {
        // Truncate the string to the maximum error length
        buf_.data()[bufErrorLength] = buf_.last() = '\0';

        FatalIOErrorInFunction(*this)
            << "problem while reading string '" << buf_.cdata() << "...' after "
            << buf_.size() << " characters\n"
            << exit(FatalIOError);

        return *this;
    }

    if (buf_.empty())
    {
        FatalIOErrorInFunction(*this)
            << "invalid first character found : " << c
            << exit(FatalIOError);
    }

    // Append end string character
    buf_.append('\0');
    str = buf_.cdata();

    // Note: check if we exited due to '}' or just !variable::valid.
    if (c != token::END_BLOCK)
    {
        putback(c);
    }

    return *this;
}


Foam::Istream& Foam::ISstream::readVerbatim(verbatimString& str)
{
    buf_.clear();

    char c;

    while (get(c))
    {
        if (c == token::HASH)
        {
            char nextC;
            get(nextC);
            if (nextC == token::END_BLOCK)
            {
                buf_.append('\0');
                str = buf_.cdata();
                return *this;
            }
            else
            {
                putback(nextC);
            }
        }

        buf_.append(c);
    }

    // Truncate the string to the maximum error length
    buf_.data()[bufErrorLength] = buf_.last() = '\0';

    FatalIOErrorInFunction(*this)
        << "problem while reading string \"" << buf_.cdata() << "...\""
        << exit(FatalIOError);

    return *this;
}


Foam::ISstream& Foam::ISstream::getLine(string& s, const bool continuation)
{
    getline(is_, s);
    setState(is_.rdstate());
    lineNumber_++;

    if (continuation && s.size())
    {
        while (s.back() == '\\')
        {
            string contLine;
            getline(is_, contLine);
            setState(is_.rdstate());
            lineNumber_++;
            s.pop_back();
            s += contLine;
        }
    }

    return *this;
}


Foam::Istream& Foam::ISstream::readDelimited
(
    string& str,
    const char begin,
    const char end
)
{
    str.clear();

    int listDepth = 0;
    char c;

    while (get(c))
    {
        str += c;

        if (c == begin)
        {
            listDepth++;
        }
        else if (c == end)
        {
            listDepth--;

            if (listDepth <= 0)
            {
                break;
            }
        }
    }

    if (bad() || listDepth != 0)
    {
        FatalIOErrorInFunction(*this)
            << "    problem while reading delimited string \n"
            << str.c_str() << endl
            << exit(FatalIOError);
    }

    return *this;
}


Foam::Istream& Foam::ISstream::readList(string& str)
{
    return readDelimited(str, token::BEGIN_LIST, token::END_LIST);
}


Foam::Istream& Foam::ISstream::readBlock(string& str)
{
    return readDelimited(str, token::BEGIN_BLOCK, token::END_BLOCK);
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


Foam::Istream& Foam::ISstream::read(longDoubleScalar& val)
{
    is_ >> val;
    setState(is_.rdstate());
    return *this;
}


Foam::Istream& Foam::ISstream::read(char* buf, std::streamsize count)
{
    if (format() != BINARY)
    {
        FatalIOErrorInFunction(*this)
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
