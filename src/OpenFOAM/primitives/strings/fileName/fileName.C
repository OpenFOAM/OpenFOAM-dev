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

#include "fileName.H"
#include "wordList.H"
#include "DynamicList.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::fileName::typeName = "fileName";
int Foam::fileName::debug(debug::debugSwitch(fileName::typeName, 0));
const Foam::fileName Foam::fileName::null;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileName::fileName(const wordList& lst)
{
    forAll(lst, elemI)
    {
        operator=((*this)/lst[elemI]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileType Foam::fileName::type
(
    const bool checkVariants,
    const bool followLink
) const
{
    return ::Foam::type(*this, checkVariants, followLink);
}


bool Foam::fileName::isName() const
{
    return find('/') == npos;
}


bool Foam::fileName::hasPath() const
{
    return find('/') != npos;
}


bool Foam::fileName::isAbsolute() const
{
    return !empty() && operator[](0) == '/';
}


Foam::fileName& Foam::fileName::toAbsolute()
{
    fileName& f = *this;

    if (!f.isAbsolute())
    {
        f = cwd()/f;
        f.clean();
    }

    return f;
}


bool Foam::fileName::clean()
{
    // The top slash - we are never allowed to go above it
    string::size_type top = this->find('/');

    // No slashes - nothing to do
    if (top == string::npos)
    {
        return false;
    }

    // Start with the '/' found:
    char prev = '/';
    string::size_type nChar  = top+1;
    string::size_type maxLen = this->size();

    for
    (
        string::size_type src = nChar;
        src < maxLen;
    )
    {
        char c = operator[](src++);

        if (prev == '/')
        {
            // Repeated '/' - skip it
            if (c == '/')
            {
                continue;
            }

            // Could be '/./' or '/../'
            if (c == '.')
            {
                // Found trailing '/.' - skip it
                if (src >= maxLen)
                {
                    continue;
                }


                // Peek at the next character
                char c1 = operator[](src);

                // Found '/./' - skip it
                if (c1 == '/')
                {
                    src++;
                    continue;
                }

                // It is '/..' or '/../'
                if (c1 == '.' && (src+1 >= maxLen || operator[](src+1) == '/'))
                {
                    string::size_type parent;

                    // Backtrack to find the parent directory
                    // Minimum of 3 characters:  '/x/../'
                    // Strip it, provided it is above the top point
                    if
                    (
                        nChar > 2
                     && (parent = this->rfind('/', nChar-2)) != string::npos
                     && parent >= top
                    )
                    {
                        nChar = parent + 1;   // Retain '/' from the parent
                        src += 2;
                        continue;
                    }

                    // Bad resolution, eg 'abc/../../'
                    // Retain the sequence, but move the top to avoid it being
                    // considered a valid parent later
                    top = nChar + 2;
                }
            }
        }
        operator[](nChar++) = prev = c;
    }

    // Remove trailing slash
    if (nChar > 1 && operator[](nChar-1) == '/')
    {
        nChar--;
    }

    this->resize(nChar);

    return (nChar != maxLen);
}


Foam::fileName Foam::fileName::clean() const
{
    fileName fName(*this);
    fName.clean();
    return fName;
}


Foam::word Foam::fileName::name() const
{
    size_type i = rfind('/');

    if (i == npos)
    {
        return *this;
    }
    else
    {
        return substr(i+1, npos);
    }
}


Foam::string Foam::fileName::caseName() const
{
    string cName = *this;

    const string caseStr(getEnv("FOAM_CASE"));

    const size_type i = find(caseStr);

    if (i == npos)
    {
        return cName;
    }
    else
    {
        return cName.replace(i, caseStr.size(), string("$FOAM_CASE"));
    }
}


Foam::word Foam::fileName::name(const bool noExt) const
{
    if (noExt)
    {
        size_type beg = rfind('/');
        if (beg == npos)
        {
            beg = 0;
        }
        else
        {
            ++beg;
        }

        size_type dot = rfind('.');
        if (dot != npos && dot <= beg)
        {
            dot = npos;
        }

        if (dot == npos)
        {
            return substr(beg, npos);
        }
        else
        {
            return substr(beg, dot - beg);
        }
    }
    else
    {
        return this->name();
    }
}


Foam::fileName Foam::fileName::path() const
{
    size_type i = rfind('/');

    if (i == npos)
    {
        return ".";
    }
    else if (i)
    {
        return substr(0, i);
    }
    else
    {
        return "/";
    }
}


Foam::fileName Foam::fileName::lessExt() const
{
    size_type i = find_last_of("./");

    if (i == npos || i == 0 || operator[](i) == '/')
    {
        return *this;
    }
    else
    {
        return substr(0, i);
    }
}


Foam::word Foam::fileName::ext() const
{
    size_type i = find_last_of("./");

    if (i == npos || i == 0 || operator[](i) == '/')
    {
        return word::null;
    }
    else
    {
        return substr(i+1, npos);
    }
}


Foam::wordList Foam::fileName::components(const char delimiter) const
{
    DynamicList<word> wrdList(20);

    size_type beg=0, end=0;

    while ((end = find(delimiter, beg)) != npos)
    {
        // Avoid empty element (caused by doubled slashes)
        if (beg < end)
        {
            wrdList.append(substr(beg, end-beg));
        }
        beg = end + 1;
    }

    // Avoid empty trailing element
    if (beg < size())
    {
        wrdList.append(substr(beg, npos));
    }

    // Transfer to wordList
    return wordList(move(wrdList));
}


Foam::word Foam::fileName::component
(
    const size_type cmpt,
    const char delimiter
) const
{
    return components(delimiter)[cmpt];
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::fileName::operator=(const fileName& str)
{
    string::operator=(str);
}


void Foam::fileName::operator=(fileName&& str)
{
    string::operator=(move(str));
}


void Foam::fileName::operator=(const word& str)
{
    string::operator=(str);
}


void Foam::fileName::operator=(const string& str)
{
    string::operator=(str);
    stripInvalid();
}


void Foam::fileName::operator=(const std::string& str)
{
    string::operator=(str);
    stripInvalid();
}


void Foam::fileName::operator=(const char* str)
{
    string::operator=(str);
    stripInvalid();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::fileName Foam::operator/(const string& a, const string& b)
{
    if (a.size())           // First string non-null
    {
        if (b.size())       // Second string non-null
        {
            return fileName(a + '/' + b);
        }
        else                // Second string null
        {
            return a;
        }
    }
    else                    // First string null
    {
        if (b.size())       // Second string non-null
        {
            return b;
        }
        else                // Second string null
        {
            return fileName();
        }
    }
}


// ************************************************************************* //
