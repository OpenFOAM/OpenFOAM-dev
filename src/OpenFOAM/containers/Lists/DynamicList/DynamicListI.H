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

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList()
:
    capacity_(0)
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
    const label nElem
)
:
    List<T>(nElem),
    capacity_(nElem)
{
    // We could also enforce SizeInc granularity when (!SizeMult || !SizeDiv)
    List<T>::size(0);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
    const label nElem,
    const T& a
)
:
    List<T>(nElem, a),
    capacity_(nElem)
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
    const DynamicList<T, SizeInc, SizeMult, SizeDiv>& lst
)
:
    List<T>(lst),
    capacity_(lst.size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
    DynamicList<T, SizeInc, SizeMult, SizeDiv>&& lst
)
:
    List<T>(move(lst)),
    capacity_(lst.capacity_)
{
    lst.capacity_ = 0;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
    const UList<T>& lst
)
:
    List<T>(lst),
    capacity_(lst.size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
    const UIndirectList<T>& lst
)
:
    List<T>(lst),
    capacity_(lst.size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
    List<T>&& lst
)
:
    List<T>(move(lst)),
    capacity_(List<T>::size())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::label Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::capacity()
const
{
    return capacity_;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::setCapacity
(
    const label nElem
)
{
    label nextFree = List<T>::size();
    capacity_ = nElem;

    if (nextFree > capacity_)
    {
        // Truncate addressed sizes too
        nextFree = capacity_;
    }

    // We could also enforce SizeInc granularity when (!SizeMult || !SizeDiv)

    List<T>::setSize(capacity_);
    List<T>::size(nextFree);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::reserve
(
    const label nElem
)
{
    // Allocate more capacity if necessary
    if (nElem > capacity_)
    {
        capacity_ = max
        (
            nElem,
            label(SizeInc + capacity_ * SizeMult / SizeDiv)
        );

        // Adjust allocated size, leave addressed size untouched
        label nextFree = List<T>::size();
        List<T>::setSize(capacity_);
        List<T>::size(nextFree);
    }
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::setSize
(
    const label nElem
)
{
    // Allocate more capacity if necessary
    if (nElem > capacity_)
    {
        capacity_ = max
        (
            nElem,
            label(SizeInc + capacity_ * SizeMult / SizeDiv)
        );

        List<T>::setSize(capacity_);
    }

    // Adjust addressed size
    List<T>::size(nElem);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::setSize
(
    const label nElem,
    const T& t
)
{
    label nextFree = List<T>::size();
    setSize(nElem);

    // Set new elements to constant value
    while (nextFree < nElem)
    {
        this->operator[](nextFree++) = t;
    }
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::resize
(
    const label nElem
)
{
    this->setSize(nElem);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::resize
(
    const label nElem,
    const T& t
)
{
    this->setSize(nElem, t);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::clear()
{
    List<T>::size(0);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::clearStorage()
{
    List<T>::clear();
    capacity_ = 0;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>&
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::shrink()
{
    label nextFree = List<T>::size();
    if (capacity_ > nextFree)
    {
        // Use the full list when resizing
        List<T>::size(capacity_);

        // The new size
        capacity_ = nextFree;
        List<T>::setSize(capacity_);
        List<T>::size(nextFree);
    }
    return *this;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::transfer
(
    DynamicList<T, SizeInc, SizeMult, SizeDiv>& lst
)
{
    // Take over storage as-is (without shrink), clear addressing for lst.
    capacity_ = lst.capacity_;
    lst.capacity_ = 0;
    List<T>::transfer(static_cast<List<T>&>(lst));
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::transfer(List<T>& lst)
{
    // Take over storage, clear addressing for lst.
    capacity_ = lst.size();
    List<T>::transfer(lst);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>&
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::append
(
    const T& t
)
{
    const label elemI = List<T>::size();
    setSize(elemI + 1);

    this->operator[](elemI) = t;
    return *this;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>&
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::append
(
    const UList<T>& lst
)
{
    if (this == &lst)
    {
        FatalErrorInFunction
            << "Attempted appending to self" << abort(FatalError);
    }

    label nextFree = List<T>::size();
    setSize(nextFree + lst.size());

    forAll(lst, elemI)
    {
        this->operator[](nextFree++) = lst[elemI];
    }
    return *this;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>&
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::append
(
    const UIndirectList<T>& lst
)
{
    label nextFree = List<T>::size();
    setSize(nextFree + lst.size());

    forAll(lst, elemI)
    {
        this->operator[](nextFree++) = lst[elemI];
    }
    return *this;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline T Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::remove()
{
    const label elemI = List<T>::size() - 1;

    if (elemI < 0)
    {
        FatalErrorInFunction
            << "List is empty" << abort(FatalError);
    }

    const T& val = List<T>::operator[](elemI);

    List<T>::size(elemI);

    return val;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline T& Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator()
(
    const label elemI
)
{
    if (elemI >= List<T>::size())
    {
        setSize(elemI + 1);
    }

    return this->operator[](elemI);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    const T& t
)
{
    UList<T>::operator=(t);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    const DynamicList<T, SizeInc, SizeMult, SizeDiv>& lst
)
{
    if (this == &lst)
    {
        FatalErrorInFunction
            << "Attempted assignment to self" << abort(FatalError);
    }

    if (capacity_ >= lst.size())
    {
        // Can copy w/o reallocating, match initial size to avoid reallocation
        List<T>::size(lst.size());
        List<T>::operator=(lst);
    }
    else
    {
        // Make everything available for the copy operation
        List<T>::size(capacity_);

        List<T>::operator=(lst);
        capacity_ = List<T>::size();
    }
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    DynamicList<T, SizeInc, SizeMult, SizeDiv>&& lst
)
{
    if (this == &lst)
    {
        FatalErrorInFunction
            << "Attempted assignment to self" << abort(FatalError);
    }

    List<T>::operator=(move(lst));
    capacity_ = lst.capacity_;
    lst.capacity_ = 0;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    List<T>&& lst
)
{
    if (this == &lst)
    {
        FatalErrorInFunction
            << "Attempted assignment to self" << abort(FatalError);
    }

    List<T>::operator=(move(lst));
    capacity_ = List<T>::size();
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    const UList<T>& lst
)
{
    if (capacity_ >= lst.size())
    {
        // Can copy w/o reallocating, match initial size to avoid reallocation
        List<T>::size(lst.size());
        List<T>::operator=(lst);
    }
    else
    {
        // Make everything available for the copy operation
        List<T>::size(capacity_);

        List<T>::operator=(lst);
        capacity_ = List<T>::size();
    }
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    const UIndirectList<T>& lst
)
{
    if (capacity_ >= lst.size())
    {
        // Can copy w/o reallocating, match initial size to avoid reallocation
        List<T>::size(lst.size());
        List<T>::operator=(lst);
    }
    else
    {
        // Make everything available for the copy operation
        List<T>::size(capacity_);

        List<T>::operator=(lst);
        capacity_ = List<T>::size();
    }
}


// * * * * * * * * * * * * * * STL Member Functions  * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
typename Foam::UList<T>::iterator
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::erase
(
    typename UList<T>::iterator curIter
)
{
    typename Foam::UList<T>::iterator iter = curIter;
    typename Foam::UList<T>::iterator nextIter = curIter;

    if (iter != this->end())
    {
        ++iter;

        while (iter != this->end())
        {
            *nextIter++ = *iter++;
        }

        this->setSize(this->size() - 1);
    }

    return curIter;
}


// ************************************************************************* //
