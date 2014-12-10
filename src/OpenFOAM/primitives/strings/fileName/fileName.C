/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
#include "debug.H"
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

Foam::fileName::Type Foam::fileName::type() const
{
    return ::Foam::type(*this);
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


//
// * remove repeated slashes
//       /abc////def        -->   /abc/def
//
// * remove '/./'
//       /abc/def/./ghi/.   -->   /abc/def/./ghi
//       abc/def/./         -->   abc/def
//
// * remove '/../'
//       /abc/def/../ghi/jkl/nmo/..   -->   /abc/ghi/jkl
//       abc/../def/ghi/../jkl        -->   abc/../def/jkl
//
// * remove trailing '/'
//
bool Foam::fileName::clean()
{
    // the top slash - we are never allowed to go above it
    register string::size_type top = this->find('/');

    // no slashes - nothing to do
    if (top == string::npos)
    {
        return false;
    }

    // start with the '/' found:
    register char prev = '/';
    register string::size_type nChar  = top+1;
    register string::size_type maxLen = this->size();

    for
    (
        register string::size_type src = nChar;
        src < maxLen;
        /*nil*/
    )
    {
        register char c = operator[](src++);

        if (prev == '/')
        {
            // repeated '/' - skip it
            if (c == '/')
            {
                continue;
            }

            // could be '/./' or '/../'
            if (c == '.')
            {
                // found trailing '/.' - skip it
                if (src >= maxLen)
                {
                    continue;
                }


                // peek at the next character
                register char c1 = operator[](src);

                // found '/./' - skip it
                if (c1 == '/')
                {
                    src++;
                    continue;
                }

                // it is '/..' or '/../'
                if (c1 == '.' && (src+1 >= maxLen || operator[](src+1) == '/'))
                {
                    string::size_type parent;

                    // backtrack to find the parent directory
                    // minimum of 3 characters:  '/x/../'
                    // strip it, provided it is above the top point
                    if
                    (
                        nChar > 2
                     && (parent = this->rfind('/', nChar-2)) != string::npos
                     && parent >= top
                    )
                    {
                        nChar = parent + 1;   // retain '/' from the parent
                        src += 2;
                        continue;
                    }

                    // bad resolution, eg 'abc/../../'
                    // retain the sequence, but move the top to avoid it being
                    // considered a valid parent later
                    top = nChar + 2;
                }
            }
        }
        operator[](nChar++) = prev = c;
    }

    // remove trailing slash
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



//  Return file name (part beyond last /)
//
//  behaviour compared to /usr/bin/basename:
//    input           name()          basename
//    -----           ------          --------
//    "foo"           "foo"           "foo"
//    "/foo"          "foo"           "foo"
//    "foo/bar"       "bar"           "bar"
//    "/foo/bar"      "bar"           "bar"
//    "/foo/bar/"     ""              "bar"
//
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


//  Return directory path name (part before last /)
//
//  behaviour compared to /usr/bin/dirname:
//    input           path()          dirname
//    -----           ------          -------
//    "foo"           "."             "."
//    "/foo"          "/"             "foo"
//    "foo/bar"       "foo"           "foo"
//    "/foo/bar"      "/foo"          "/foo"
//    "/foo/bar/"     "/foo/bar/"     "/foo"
//
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


//  Return file name without extension (part before last .)
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


//  Return file name extension (part after last .)
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


// Return the components of the file name as a wordList
// note that concatenating the components will not necessarily retrieve
// the original input fileName
//
//  behaviour
//    input           components()
//    -----           ------
//    "foo"           1("foo")
//    "/foo"          1("foo")
//    "foo/bar"       2("foo", "bar")
//    "/foo/bar"      2("foo", "bar")
//    "/foo/bar/"     2("foo", "bar")
//
Foam::wordList Foam::fileName::components(const char delimiter) const
{
    DynamicList<word> wrdList(20);

    size_type beg=0, end=0;

    while ((end = find(delimiter, beg)) != npos)
    {
        // avoid empty element (caused by doubled slashes)
        if (beg < end)
        {
            wrdList.append(substr(beg, end-beg));
        }
        beg = end + 1;
    }

    // avoid empty trailing element
    if (beg < size())
    {
        wrdList.append(substr(beg, npos));
    }

    // transfer to wordList
    return wordList(wrdList.xfer());
}


// Return a component of the file name
Foam::word Foam::fileName::component
(
    const size_type cmpt,
    const char delimiter
) const
{
    return components(delimiter)[cmpt];
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

const Foam::fileName& Foam::fileName::operator=(const fileName& str)
{
    string::operator=(str);
    return *this;
}


const Foam::fileName& Foam::fileName::operator=(const word& str)
{
    string::operator=(str);
    return *this;
}


const Foam::fileName& Foam::fileName::operator=(const string& str)
{
    string::operator=(str);
    stripInvalid();
    return *this;
}


const Foam::fileName& Foam::fileName::operator=(const std::string& str)
{
    string::operator=(str);
    stripInvalid();
    return *this;
}


const Foam::fileName& Foam::fileName::operator=(const char* str)
{
    string::operator=(str);
    stripInvalid();
    return *this;
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
