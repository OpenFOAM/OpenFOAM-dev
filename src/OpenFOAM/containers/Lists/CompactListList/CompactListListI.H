/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "CompactListList.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::CompactListList<T>::CompactListList()
:
    UCompactListList<T>(),
    offsets_(label(1), label(0)),
    m_()
{
    UCompactListList<T>::shallowCopy(UCompactListList<T>(offsets_, m_));
}


template<class T>
inline Foam::CompactListList<T>::CompactListList
(
    const UList<label>& offsets,
    const UList<T>& m
)
:
    UCompactListList<T>(offsets, m),
    offsets_(offsets),
    m_(m)
{}


template<class T>
Foam::CompactListList<T>::CompactListList(CompactListList<T>& lst, bool reuse)
:
    UCompactListList<T>(),
    offsets_(lst.offsets_, reuse),
    m_(lst.m_, reuse)
{
    UCompactListList<T>::shallowCopy(UCompactListList<T>(offsets_, m_));
}


template<class T>
Foam::CompactListList<T>::CompactListList(CompactListList<T>&& lst)
{
    transfer(lst);
}


template<class T>
inline Foam::CompactListList<T>::CompactListList
(
    const label mRows,
    const label nData,
    const T& t
)
:
    UCompactListList<T>(),
    offsets_(mRows + 1, 0),
    m_(nData, t)
{
    UCompactListList<T>::shallowCopy(UCompactListList<T>(offsets_, m_));
}


template<class T>
inline Foam::autoPtr<Foam::CompactListList<T>>
Foam::CompactListList<T>::clone() const
{
    return autoPtr<CompactListList<T>>
    (
        new CompactListList<T>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const Foam::CompactListList<T>&
Foam::CompactListList<T>::null()
{
    return NullObjectRef<CompactListList<T>>();
}


template<class T>
inline void Foam::CompactListList<T>::resize(const label mRows)
{
    this->setSize(mRows);
}


template<class T>
inline void Foam::CompactListList<T>::resize
(
    const label mRows,
    const label nData
)
{
    this->setSize(mRows, nData);
}


template<class T>
inline void Foam::CompactListList<T>::resize
(
    const label mRows,
    const label nData,
    const T& t
)
{
    this->setSize(mRows, nData, t);
}


template<class T>
inline void Foam::CompactListList<T>::resize
(
    const labelUList& rowSizes
)
{
    this->setSize(rowSizes);
}


template<class T>
template<class T2>
void Foam::CompactListList<T>::resize(const UList<T2>& ll)
{
    this->setSize(ll);
}


template<class T>
template<class T2>
void Foam::CompactListList<T>::resize(const UIndirectList<T2>& ll)
{
    this->setSize(ll);
}


// ************************************************************************* //
