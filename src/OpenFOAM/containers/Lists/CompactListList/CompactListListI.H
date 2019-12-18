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

#include "ListOps.H"
#include "SubList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Container>
inline Foam::CompactListList<T, Container>::CompactListList()
:
    size_(0)
{}


template<class T, class Container>
inline Foam::CompactListList<T, Container>::CompactListList
(
    const label mRows,
    const label nData
)
:
    size_(mRows),
    offsets_(mRows+1, 0),
    m_(nData)
{}


template<class T, class Container>
inline Foam::CompactListList<T, Container>::CompactListList
(
    const label mRows,
    const label nData,
    const T& t
)
:
    size_(mRows),
    offsets_(mRows+1, 0),
    m_(nData, t)
{}


template<class T, class Container>
inline Foam::autoPtr<Foam::CompactListList<T, Container>>
Foam::CompactListList<T, Container>::clone() const
{
    return autoPtr<CompactListList<T, Container>>
    (
        new CompactListList<T, Container>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class Container>
inline const Foam::CompactListList<T, Container>&
Foam::CompactListList<T, Container>::null()
{
    return NullObjectRef<CompactListList<T, Container>>();
}


template<class T, class Container>
inline Foam::label Foam::CompactListList<T, Container>::size() const
{
    return size_;
}


template<class T, class Container>
inline bool Foam::CompactListList<T, Container>::empty() const
{
    return !size_;
}


template<class T, class Container>
inline const Foam::List<Foam::label>&
Foam::CompactListList<T, Container>::offsets() const
{
    return offsets_;
}


template<class T, class Container>
inline Foam::List<Foam::label>& Foam::CompactListList<T, Container>::offsets()
{
    return offsets_;
}


template<class T, class Container>
inline const Foam::List<T>& Foam::CompactListList<T, Container>::m()
const
{
    return m_;
}


template<class T, class Container>
inline Foam::List<T>& Foam::CompactListList<T, Container>::m()
{
    return m_;
}


template<class T, class Container>
inline Foam::label Foam::CompactListList<T, Container>::index
(
    const label i,
    const label j
) const
{
    return offsets_[i] + j;
}


template<class T, class Container>
inline Foam::label Foam::CompactListList<T, Container>::whichRow(const label i)
const
{
    if (i < 0 || i >= m_.size())
    {
        FatalErrorInFunction
            << "Index " << i << " outside 0.." << m_.size()
            << abort(FatalError);
    }

    return findLower(offsets_, i+1);
}


template<class T, class Container>
inline Foam::label Foam::CompactListList<T, Container>::whichColumn
(
    const label row,
    const label i
) const
{
    return i - index(row, 0);
}


template<class T, class Container>
inline void Foam::CompactListList<T, Container>::resize(const label mRows)
{
    this->setSize(mRows);
}


template<class T, class Container>
inline void Foam::CompactListList<T, Container>::resize
(
    const label mRows,
    const label nData
)
{
    this->setSize(mRows, nData);
}


template<class T, class Container>
inline void Foam::CompactListList<T, Container>::resize
(
    const label mRows,
    const label nData,
    const T& t
)
{
    this->setSize(mRows, nData, t);
}


template<class T, class Container>
inline void Foam::CompactListList<T, Container>::resize
(
    const labelUList& rowSizes
)
{
    this->setSize(rowSizes);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Container>
inline Foam::UList<T> Foam::CompactListList<T, Container>::operator[]
(
    const label i
)
{
    label start = offsets_[i];
    return UList<T>(m_.begin() + start, offsets_[i+1] - start);
}


template<class T, class Container>
inline const Foam::UList<T>
Foam::CompactListList<T, Container>::operator[]
(
    const label i
) const
{
    label start = offsets_[i];
    return UList<T>
    (
        const_cast<T*>(m_.begin() + start),
        offsets_[i+1] - start
    );
}


template<class T, class Container>
inline T& Foam::CompactListList<T, Container>::operator()
(
    const label i,
    const label j
)
{
    return m_[index(i, j)];
}


template<class T, class Container>
inline const T& Foam::CompactListList<T, Container>::operator()
(
    const label i,
    const label j
) const
{
    return m_[index(i, j)];
}


template<class T, class Container>
inline void Foam::CompactListList<T, Container>::operator=
(
    CompactListList<T, Container>&& cll
)
{
    transfer(cll);
}


template<class T, class Container>
inline void Foam::CompactListList<T, Container>::operator=(const T& t)
{
    m_ = t;
}


// ************************************************************************* //
