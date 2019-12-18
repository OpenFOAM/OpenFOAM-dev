/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::PackedBoolList::PackedBoolList()
:
    PackedList<1>()
{}


inline Foam::PackedBoolList::PackedBoolList(const label size)
:
    PackedList<1>(size)
{}


inline Foam::PackedBoolList::PackedBoolList
(
    const label size,
    const bool val
)
:
    PackedList<1>(size, (val ? 1u : 0u))
{}


inline Foam::PackedBoolList::PackedBoolList(const PackedBoolList& lst)
:
    PackedList<1>(lst)
{}


inline Foam::PackedBoolList::PackedBoolList(const PackedList<1>& lst)
:
    PackedList<1>(lst)
{}


inline Foam::PackedBoolList::PackedBoolList(const Foam::UList<bool>& lst)
:
    PackedList<1>()
{
    operator=(lst);
}


inline Foam::PackedBoolList::PackedBoolList(const labelUList& indices)
:
    PackedList<1>(indices.size(), 0u)
{
    set(indices);
}


inline Foam::PackedBoolList::PackedBoolList(const UIndirectList<label>& indices)
:
    PackedList<1>(indices.size(), 0u)
{
    set(indices);
}


inline Foam::autoPtr<Foam::PackedBoolList>
Foam::PackedBoolList::clone() const
{
    return autoPtr<PackedBoolList>(new PackedBoolList(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::PackedBoolList::transfer(PackedBoolList& lst)
{
    PackedList<1>::transfer(static_cast<PackedList<1>&>(lst));
}


inline void Foam::PackedBoolList::transfer(PackedList<1>& lst)
{
    PackedList<1>::transfer(lst);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::PackedBoolList::operator=(const bool val)
{
    PackedList<1>::operator=(val);
}


inline void Foam::PackedBoolList::operator=(const PackedBoolList& lst)
{
    PackedList<1>::operator=(lst);
}


inline void Foam::PackedBoolList::operator=(const PackedList<1>& lst)
{
    PackedList<1>::operator=(lst);
}


inline void Foam::PackedBoolList::operator=(const labelUList& indices)
{
    clear();
    set(indices);
}


inline void Foam::PackedBoolList::operator=(const UIndirectList<label>& indices)
{
    clear();
    set(indices);
}


inline Foam::PackedBoolList
Foam::PackedBoolList::operator~() const
{
    PackedBoolList result(*this);
    result.flip();

    return result;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator&=(const PackedList<1>& lst)
{
    subset(lst);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator&=(const labelUList& indices)
{
    subset(indices);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator&=(const UIndirectList<label>& indices)
{
    subset(indices);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator|=(const PackedList<1>& lst)
{
    set(lst);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator|=(const labelUList& indices)
{
    set(indices);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator|=(const UIndirectList<label>& indices)
{
    set(indices);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator+=(const PackedList<1>& lst)
{
    return operator|=(lst);
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator+=(const labelUList& indices)
{
    return operator|=(indices);
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator+=(const UIndirectList<label>& indices)
{
    return operator|=(indices);
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator-=(const PackedList<1>& lst)
{
    unset(lst);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator-=(const labelUList& indices)
{
    unset(indices);
    return *this;
}


inline Foam::PackedBoolList&
Foam::PackedBoolList::operator-=(const UIndirectList<label>& indices)
{
    unset(indices);
    return *this;
}


// ************************************************************************* //
