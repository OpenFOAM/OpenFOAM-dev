/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include <iostream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::string::string()
{}


inline Foam::string::string(const std::string& str)
:
    std::string(str)
{}


// Copy character array
inline Foam::string::string(const char* str)
:
    std::string(str)
{}


// Construct from a given number of characters in a character array
inline Foam::string::string(const char* str, const size_type len)
:
    std::string(str, len)
{}


// Construct from a single character
inline Foam::string::string(const char c)
:
    std::string(1, c)
{}


inline Foam::string::string(const size_type len, const char c)
:
    std::string(len, c)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class String>
inline bool Foam::string::valid(const string& str)
{
    for (const_iterator iter = str.begin(); iter != str.end(); ++iter)
    {
        if (!String::valid(*iter))
        {
            return false;
        }
    }
    return true;
}


template<class String>
inline bool Foam::string::stripInvalid(string& str)
{
    if (!valid<String>(str))
    {
        size_type nValid = 0;
        iterator iter2 = str.begin();

        for
        (
            const_iterator iter1 = iter2;
            iter1 != const_cast<const string&>(str).end();
            iter1++
        )
        {
            char c = *iter1;

            if (String::valid(c))
            {
                *iter2 = c;
                ++iter2;
                ++nValid;
            }
        }

        str.resize(nValid);

        return true;
    }

    return false;
}


template<class String>
inline bool Foam::string::meta(const string& str, const char quote)
{
    int escaped = 0;
    for (const_iterator iter = str.begin(); iter != str.end(); ++iter)
    {
        if (quote && *iter == quote)
        {
            escaped ^= 1;  // toggle state
        }
        else if (escaped)
        {
            escaped = false;
        }
        else if (String::meta(*iter))
        {
            return true;
        }
    }
    return false;
}


template<class String>
inline Foam::string
Foam::string::quotemeta(const string& str, const char quote)
{
    if (!quote)
    {
        return str;
    }

    string sQuoted;
    sQuoted.reserve(2*str.length());

    int escaped = 0;
    for (const_iterator iter = str.begin(); iter != str.end(); ++iter)
    {
        if (*iter == quote)
        {
            escaped ^= 1;  // toggle state
        }
        else if (escaped)
        {
            escaped = 0;
        }
        else if (String::meta(*iter))
        {
            sQuoted += quote;
        }

        sQuoted += *iter;
    }

    sQuoted.resize(sQuoted.length());

    return sQuoted;
}


template<class String>
inline String Foam::string::validate(const string& str)
{
    string ss = str;
    stripInvalid<String>(ss);
    return ss;
}

inline bool Foam::string::match(const std::string& str) const
{
    // check as string
    return (str == *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::string Foam::string::operator()
(
    const size_type i,
    const size_type n
) const
{
    return substr(i, n);
}


inline Foam::string Foam::string::operator()(const size_type n) const
{
    return substr(0, n);
}


inline unsigned Foam::string::hash::operator()
(
    const string& key,
    unsigned seed
) const
{
    return Hasher(key.data(), key.size(), seed);
}

// ************************************************************************* //
