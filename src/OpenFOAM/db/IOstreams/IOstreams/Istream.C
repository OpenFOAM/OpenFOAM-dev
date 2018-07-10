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

#include "Istream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Istream::putBack(const token& t)
{
    if (bad())
    {
        FatalIOErrorInFunction(*this)
            << "Attempt to put back onto bad stream"
            << exit(FatalIOError);
    }
    else if (putBack_)
    {
        FatalIOErrorInFunction(*this)
            << "Attempt to put back another token"
            << exit(FatalIOError);
    }
    else
    {
        putBackToken_ = t;
        putBack_ = true;
    }
}


bool Foam::Istream::getBack(token& t)
{
    if (bad())
    {
        FatalIOErrorInFunction(*this)
            << "Attempt to get back from bad stream"
            << exit(FatalIOError);
    }
    else if (putBack_)
    {
        t = putBackToken_;
        putBack_ = false;
        return true;
    }

    return false;
}


bool Foam::Istream::peekBack(token& t)
{
    if (putBack_)
    {
        t = putBackToken_;
    }
    else
    {
        t = token::undefinedToken;
    }

    return putBack_;
}


Foam::Istream& Foam::Istream::readBegin(const char* funcName)
{
    token delimiter(*this);
    if (delimiter != token::BEGIN_LIST)
    {
        setBad();
        FatalIOErrorInFunction(*this)
            << "Expected a '" << token::BEGIN_LIST
            << "' while reading " << funcName
            << ", found " << delimiter.info()
            << exit(FatalIOError);
    }

    return *this;
}


Foam::Istream& Foam::Istream::readEnd(const char* funcName)
{
    token delimiter(*this);
    if (delimiter != token::END_LIST)
    {
        setBad();
        FatalIOErrorInFunction(*this)
            << "Expected a '" << token::END_LIST
            << "' while reading " << funcName
            << ", found " << delimiter.info()
            << exit(FatalIOError);
    }

    return *this;
}


Foam::Istream& Foam::Istream::readEndBegin(const char* funcName)
{
    readEnd(funcName);
    return readBegin(funcName);
}


char Foam::Istream::readBeginList(const char* funcName)
{
    token delimiter(*this);

    if (delimiter != token::BEGIN_LIST && delimiter != token::BEGIN_BLOCK)
    {
        setBad();
        FatalIOErrorInFunction(*this)
            << "Expected a '" << token::BEGIN_LIST
            << "' or a '" << token::BEGIN_BLOCK
            << "' while reading " << funcName
            << ", found " << delimiter.info()
            << exit(FatalIOError);

        return '\0';
    }

    return delimiter.pToken();
}


char Foam::Istream::readEndList(const char* funcName)
{
    token delimiter(*this);

    if (delimiter != token::END_LIST && delimiter != token::END_BLOCK)
    {
        setBad();
        FatalIOErrorInFunction(*this)
            << "Expected a '" << token::END_LIST
            << "' or a '" << token::END_BLOCK
            << "' while reading " << funcName
            << ", found " << delimiter.info()
            << exit(FatalIOError);

        return '\0';
    }

    return delimiter.pToken();
}


Foam::Istream& Foam::Istream::operator()() const
{
    if (!good())
    {
        check("Istream::operator()");
        FatalIOError.exit();
    }

    return const_cast<Istream&>(*this);
}


// ************************************************************************* //
