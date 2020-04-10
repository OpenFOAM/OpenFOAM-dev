/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "wordAndDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wordAndDictionary::wordAndDictionary()
:
    Tuple2<word, dictionary>()
{}


Foam::wordAndDictionary::wordAndDictionary(Istream& is)
:
    Tuple2<word, dictionary>()
{
    is >> *this;

    is.check("wordAndDictionary::wordAndDictionary(Istream& is)");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wordAndDictionary::~wordAndDictionary()
{}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, wordAndDictionary& wd)
{
    wd.first() = word(is);
    wd.second().clear();

    token t(is);
    is.putBack(t);

    if (t.isPunctuation() && t.pToken() == token::BEGIN_BLOCK)
    {
        dictionary d(is);
        wd.second().transfer(d);
    }

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const wordAndDictionary& wd)
{
    os << wd.first();

    if (!wd.second().empty())
    {
        os << token::SPACE << wd.second();
    }

    return os;
}


// ************************************************************************* //
