/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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


/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

const char* const Foam::string::typeName = "string";
int Foam::string::debug(Foam::debug::debugSwitch(string::typeName, 0));
const Foam::string Foam::string::null;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Count and return the number of a given character in the string
Foam::string::size_type Foam::string::count(const char c) const
{
    register size_type cCount = 0;

    for (const_iterator iter = begin(); iter != end(); ++iter)
    {
        if (*iter == c)
        {
            ++cCount;
        }
    }

    return cCount;
}


// Replace first occurence of sub-string oldStr with newStr
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


// Replace all occurences of sub-string oldStr with newStr
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


// Remove repeated characters returning true if string changed
bool Foam::string::removeRepeated(const char character)
{
    bool changed = false;

    if (character && find(character) != npos)
    {
        register string::size_type nChar=0;
        iterator iter2 = begin();

        register char prev = 0;

        for
        (
            string::const_iterator iter1 = iter2;
            iter1 != end();
            iter1++
        )
        {
            register char c = *iter1;

            if (prev == c && c == character)
            {
                changed = true;
            }
            else
            {
                *iter2 = prev = c;
                ++iter2;
                ++nChar;
            }
        }
        resize(nChar);
    }

    return changed;
}


// Return string with repeated characters removed
Foam::string Foam::string::removeRepeated(const char character) const
{
    string str(*this);
    str.removeRepeated(character);
    return str;
}


// Remove trailing character returning true if string changed
bool Foam::string::removeTrailing(const char character)
{
    bool changed = false;

    string::size_type nChar = size();
    if (character && nChar > 1 && operator[](nChar-1) == character)
    {
        resize(nChar-1);
        changed = true;
    }

    return changed;
}


// Return string with trailing character removed
Foam::string Foam::string::removeTrailing(const char character) const
{
    string str(*this);
    str.removeTrailing(character);
    return str;
}


// ************************************************************************* //
