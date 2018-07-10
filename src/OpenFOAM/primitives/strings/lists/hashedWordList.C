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

#include "hashedWordList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::hashedWordList::rehash()
{
    indices_.clear();
    forAll(*this, i)
    {
        indices_.insert(List<word>::operator[](i), i);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hashedWordList::hashedWordList()
:
    List<word>()
{}


Foam::hashedWordList::hashedWordList(const UList<word>& names)
:
    List<word>(names)
{
    rehash();
}


Foam::hashedWordList::hashedWordList(const hashedWordList& names)
:
    List<word>(static_cast<const UList<word>&>(names))
{
    rehash();
}


Foam::hashedWordList::hashedWordList(const Xfer<List<word>>& names)
:
    List<word>(names)
{
    rehash();
}


Foam::hashedWordList::hashedWordList
(
    const label nNames,
    const char** names
)
:
    List<word>(nNames)
{
    forAll(*this, i)
    {
        List<word>::operator[](i) = names[i];
    }

    rehash();
}


Foam::hashedWordList::hashedWordList
(
    const char** names
)
{
    // count names
    label nNames = 0;
    for (unsigned i = 0; names[i] && *(names[i]); ++i)
    {
        ++nNames;
    }

    List<word>::setSize(nNames);
    forAll(*this, i)
    {
        List<word>::operator[](i) = names[i];
    }

    rehash();
}


Foam::hashedWordList::hashedWordList(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::hashedWordList::clear()
{
    List<word>::clear();
    indices_.clear();
}


void Foam::hashedWordList::append(const word& name)
{
    const label idx = size();
    List<word>::append(name);
    indices_.insert(name, idx);
}


void Foam::hashedWordList::transfer(List<word>& lst)
{
    List<word>::transfer(lst);
    rehash();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, hashedWordList& lst)
{
    is  >> static_cast<List<word>&>(lst);
    lst.rehash();

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const hashedWordList& lst)
{
    os  << static_cast<const List<word>&>(lst);
    return os;
}


// ************************************************************************* //
