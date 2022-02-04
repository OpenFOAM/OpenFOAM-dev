/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "UCompactListList.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::UCompactListList<T>::UCompactListList()
:
    offsets_(nullptr, 1),
    m_()
{}


template<class T>
inline Foam::UCompactListList<T>::UCompactListList
(
    const UList<label>& offsets,
    const UList<T>& m
)
:
    offsets_(offsets),
    m_(m)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const Foam::UCompactListList<T>&
Foam::UCompactListList<T>::null()
{
    return NullObjectRef<UCompactListList<T>>();
}


template<class T>
inline Foam::label Foam::UCompactListList<T>::size() const
{
    return offsets_.size() - 1;
}


template<class T>
inline bool Foam::UCompactListList<T>::empty() const
{
    return !size();
}


template<class T>
inline const Foam::UList<Foam::label>&
Foam::UCompactListList<T>::offsets() const
{
    return offsets_;
}


template<class T>
inline Foam::UList<Foam::label>& Foam::UCompactListList<T>::offsets()
{
    return offsets_;
}


template<class T>
inline const Foam::UList<T>& Foam::UCompactListList<T>::m()
const
{
    return m_;
}


template<class T>
inline Foam::UList<T>& Foam::UCompactListList<T>::m()
{
    return m_;
}


template<class T>
inline void Foam::UCompactListList<T>::shallowCopy
(
    const UCompactListList<T>& l
)
{
    offsets_.shallowCopy(l.offsets_);
    m_.shallowCopy(l.m_);
}


template<class T>
inline void Foam::UCompactListList<T>::deepCopy
(
    const UCompactListList<T>& l
)
{
    offsets_.deepCopy(l.offsets_);
    m_.deepCopy(l.m_);
}


template<class T>
inline Foam::label Foam::UCompactListList<T>::index
(
    const label i,
    const label j
) const
{
    return offsets_[i] + j;
}


template<class T>
inline Foam::label Foam::UCompactListList<T>::whichRow(const label i) const
{
    if (i < 0 || i >= m_.size())
    {
        FatalErrorInFunction
            << "Index " << i << " outside 0.." << m_.size()
            << abort(FatalError);
    }

    return findLower(offsets_, i+1);
}


template<class T>
inline Foam::label Foam::UCompactListList<T>::whichColumn
(
    const label row,
    const label i
) const
{
    return i - index(row, 0);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline Foam::UList<T> Foam::UCompactListList<T>::operator[]
(
    const label i
)
{
    label start = offsets_[i];
    return UList<T>(m_.begin() + start, offsets_[i+1] - start);
}


template<class T>
inline const Foam::UList<T>
Foam::UCompactListList<T>::operator[]
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


template<class T>
inline T& Foam::UCompactListList<T>::operator()
(
    const label i,
    const label j
)
{
    return m_[index(i, j)];
}


template<class T>
inline const T& Foam::UCompactListList<T>::operator()
(
    const label i,
    const label j
) const
{
    return m_[index(i, j)];
}


template<class T>
inline void Foam::UCompactListList<T>::operator=(const T& t)
{
    m_ = t;
}


// ************************************************************************* //
