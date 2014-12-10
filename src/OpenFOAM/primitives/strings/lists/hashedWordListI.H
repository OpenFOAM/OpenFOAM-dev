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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::hashedWordList::contains(const word& name) const
{
    return indices_.found(name);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::hashedWordList::operator=(const UList<word>& lst)
{
    List<word>::operator=(lst);
    rehash();
}


inline void Foam::hashedWordList::operator=(const hashedWordList& lst)
{
    operator=(static_cast<const UList<word>&>(lst));
}


inline const Foam::word& Foam::hashedWordList::operator[]
(
    const label index
) const
{
    return List<word>::operator[](index);
}


// could return -1 instead of bombing out
inline Foam::label Foam::hashedWordList::operator[](const word& name) const
{
    return indices_[name];
}


// ************************************************************************* //
