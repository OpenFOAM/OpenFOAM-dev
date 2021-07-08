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

#include "string.H"
#include "stringOps.H"
#include "UList.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

const char* const Foam::string::typeName = "string";
int Foam::string::debug(Foam::debug::debugSwitch(string::typeName, 0));
const Foam::string Foam::string::null;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::string::string(const UList<char>& str)
:
    std::string(str.begin(), str.end())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::string::size_type Foam::string::count(const char c) const
{
    size_type cCount = 0;

    for (const_iterator iter = begin(); iter != end(); ++iter)
    {
        if (*iter == c)
        {
            ++cCount;
        }
    }

    return cCount;
}


Foam::string& Foam::string::replace
(
    const string& oldStr,
    const string& newStr,
    size_type start
)
{
    size_type newStart = start;

    if ((newStart = find(oldStr, newStart)) != npos)
    {
        std::string::replace(newStart, oldStr.size(), newStr);
    }

    return *this;
}


Foam::string& Foam::string::replaceAll
(
    const string& oldStr,
    const string& newStr,
    size_type start
)
{
    if (oldStr.size())
    {
        size_type newStart = start;

        while ((newStart = find(oldStr, newStart)) != npos)
        {
            std::string::replace(newStart, oldStr.size(), newStr);
            newStart += newStr.size();
        }
    }

    return *this;
}


Foam::string& Foam::string::expand(const bool allowEmpty)
{
    stringOps::inplaceExpand(*this, allowEmpty);
    return *this;
}


bool Foam::string::removeRepeated(const char character)
{
    bool changed = false;

    string::size_type n = 0;
    iterator iter2 = begin();

    char cPrev = operator[](0) + 1;

    for
    (
        string::const_iterator iter1 = iter2;
        iter1 != end();
        ++ iter1
    )
    {
        char c = *iter1;

        if (c == cPrev && c == character)
        {
            changed = true;
        }
        else
        {
            *iter2 = cPrev = c;
            ++ iter2;
            ++ n;
        }
    }

    resize(n);

    return changed;
}


Foam::string Foam::string::removeRepeated(const char character) const
{
    string str(*this);
    str.removeRepeated(character);
    return str;
}


bool Foam::string::removeTrailing(const char character)
{
    bool changed = false;

    string::size_type n = size();
    if (n >= 1 && operator[](n - 1) == character)
    {
        resize(n - 1);
        changed = true;
    }

    return changed;
}


Foam::string Foam::string::removeTrailing(const char character) const
{
    string result(*this);
    result.removeTrailing(character);
    return result;
}


bool Foam::string::removeTrailing(const string& str)
{
    bool changed = false;

    string::size_type n = size(), nStr = str.size();
    if (n >= str.size() && operator()(n - nStr, nStr) == str)
    {
        resize(n - nStr);
        changed = true;
    }

    return changed;
}


Foam::string Foam::string::removeTrailing(const string& str) const
{
    string result(*this);
    result.removeTrailing(str);
    return result;
}


void Foam::string::strip(const string& str)
{
    // Find the first character to keep
    string::size_type i0 = 0;
    while (i0 < size() && str.count(operator[](i0)) > 0)
    {
        ++ i0;
    }

    // Find one past the last character to keep
    string::size_type i1 = size();
    while (i1 > i0 && str.count(operator[](i1 - 1)) > 0)
    {
        -- i1;
    }

    // Remove leading characters by shuffling the string up
    if (i0 != 0)
    {
        for (string::size_type i = 0; i < size() - i0; ++ i)
        {
            operator[](i) = operator[](i + i0);
        }
    }

    // If removing any characters then resize the string
    if (i0 != 0 || i1 != size())
    {
        resize(i1 - i0);
    }
}


// ************************************************************************* //
