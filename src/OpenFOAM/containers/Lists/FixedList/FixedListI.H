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

#include "UList.H"
#include "SLList.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, unsigned Size>
inline Foam::FixedList<T, Size>::FixedList()
{}


template<class T, unsigned Size>
inline Foam::FixedList<T, Size>::FixedList(const T v[Size])
{
    for (register unsigned i=0; i<Size; i++)
    {
        v_[i] = v[i];
    }
}


template<class T, unsigned Size>
inline Foam::FixedList<T, Size>::FixedList(const T& t)
{
    for (register unsigned i=0; i<Size; i++)
    {
        v_[i] = t;
    }
}


template<class T, unsigned Size>
inline Foam::FixedList<T, Size>::FixedList(const UList<T>& lst)
{
    checkSize(lst.size());

    for (register unsigned i=0; i<Size; i++)
    {
        v_[i] = lst[i];
    }
}


template<class T, unsigned Size>
inline Foam::FixedList<T, Size>::FixedList(const SLList<T>& lst)
{
    checkSize(lst.size());

    register label i = 0;
    for
    (
        typename SLList<T>::const_iterator iter = lst.begin();
        iter != lst.end();
        ++iter
    )
    {
        operator[](i++) = iter();
    }
}


template<class T, unsigned Size>
inline Foam::FixedList<T, Size>::FixedList(const FixedList<T, Size>& lst)
{
    for (register unsigned i=0; i<Size; i++)
    {
        v_[i] = lst[i];
    }
}


template<class T, unsigned Size>
inline Foam::autoPtr< Foam::FixedList<T, Size> >
Foam::FixedList<T, Size>::clone() const
{
    return autoPtr< FixedList<T, Size> >(new FixedList<T, Size>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, unsigned Size>
inline const Foam::FixedList<T, Size>& Foam::FixedList<T, Size>::null()
{
    return *reinterpret_cast< FixedList<T, Size>* >(0);
}


template<class T, unsigned Size>
inline Foam::label Foam::FixedList<T, Size>::fcIndex(const label i) const
{
    return (i == Size-1 ? 0 : i+1);
}


template<class T, unsigned Size>
inline Foam::label Foam::FixedList<T, Size>::rcIndex(const label i) const
{
    return (i ? i-1 : Size-1);
}


// Check start is within valid range (0 ... size-1).
template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::checkStart(const label start) const
{
    if (start < 0 || (start && unsigned(start) >= Size))
    {
        FatalErrorIn("FixedList<T, Size>::checkStart(const label)")
            << "start " << start << " out of range 0 ... " << (Size-1)
            << abort(FatalError);
    }
}


// Check size is within valid range (0 ... size).
template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::checkSize(const label size) const
{
    if (size < 0 || unsigned(size) > Size)
    {
        FatalErrorIn("FixedList<T, Size>::checkSize(const label)")
            << "size " << size << " out of range 0 ... " << (Size)
            << abort(FatalError);
    }
}


// Check index i is within valid range (0 ... size-1)
// The check for zero-sized list is already done in static assert
template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::checkIndex(const label i) const
{
    if (i < 0 || unsigned(i) >= Size)
    {
        FatalErrorIn("FixedList<T, Size>::checkIndex(const label)")
            << "index " << i << " out of range 0 ... " << (Size-1)
            << abort(FatalError);
    }
}


template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::resize(const label s)
{
#   ifdef FULLDEBUG
    checkSize(s);
#   endif
}

template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::setSize(const label s)
{
#   ifdef FULLDEBUG
    checkSize(s);
#   endif
}

template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::transfer(const FixedList<T, Size>& lst)
{
    for (register unsigned i=0; i<Size; i++)
    {
        v_[i] = lst[i];
    }
}


template<class T, unsigned Size>
inline const T*
Foam::FixedList<T, Size>::cdata() const
{
    return v_;
}


template<class T, unsigned Size>
inline T*
Foam::FixedList<T, Size>::data()
{
    return v_;
}


template<class T, unsigned Size>
inline T& Foam::FixedList<T, Size>::first()
{
    return v_[0];
}


template<class T, unsigned Size>
inline const T& Foam::FixedList<T, Size>::first() const
{
    return v_[0];
}


template<class T, unsigned Size>
inline T& Foam::FixedList<T, Size>::last()
{
    return v_[Size-1];
}


template<class T, unsigned Size>
inline const T& Foam::FixedList<T, Size>::last() const
{
    return v_[Size-1];
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// element access
template<class T, unsigned Size>
inline T& Foam::FixedList<T, Size>::operator[](const label i)
{
#   ifdef FULLDEBUG
    checkIndex(i);
#   endif
    return v_[i];
}


// const element access
template<class T, unsigned Size>
inline const T& Foam::FixedList<T, Size>::operator[](const label i) const
{
#   ifdef FULLDEBUG
    checkIndex(i);
#   endif
    return v_[i];
}


template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::operator=(const T lst[Size])
{
    for (register unsigned i=0; i<Size; i++)
    {
        v_[i] = lst[i];
    }
}

template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::operator=(const UList<T>& lst)
{
    checkSize(lst.size());

    for (register unsigned i=0; i<Size; i++)
    {
        v_[i] = lst[i];
    }
}

template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::operator=(const SLList<T>& lst)
{
    checkSize(lst.size());

    register label i = 0;
    for
    (
        typename SLList<T>::const_iterator iter = lst.begin();
        iter != lst.end();
        ++iter
    )
    {
        operator[](i++) = iter();
    }
}

template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::operator=(const T& t)
{
    for (register unsigned i=0; i<Size; i++)
    {
        v_[i] = t;
    }
}


// * * * * * * * * * * * * * * STL Member Functions  * * * * * * * * * * * * //

template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::iterator
Foam::FixedList<T, Size>::begin()
{
    return v_;
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::const_iterator
Foam::FixedList<T, Size>::begin() const
{
    return v_;
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::const_iterator
Foam::FixedList<T, Size>::cbegin() const
{
    return v_;
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::iterator
Foam::FixedList<T, Size>::end()
{
    return &v_[Size];
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::const_iterator
Foam::FixedList<T, Size>::end() const
{
    return &v_[Size];
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::const_iterator
Foam::FixedList<T, Size>::cend() const
{
    return &v_[Size];
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::iterator
Foam::FixedList<T, Size>::rbegin()
{
    return &v_[Size-1];
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::const_iterator
Foam::FixedList<T, Size>::rbegin() const
{
    return &v_[Size-1];
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::const_iterator
Foam::FixedList<T, Size>::crbegin() const
{
    return &v_[Size-1];
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::iterator
Foam::FixedList<T, Size>::rend()
{
    return &v_[-1];
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::const_iterator
Foam::FixedList<T, Size>::rend() const
{
    return &v_[-1];
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::const_iterator
Foam::FixedList<T, Size>::crend() const
{
    return &v_[-1];
}


template<class T, unsigned Size>
inline Foam::label Foam::FixedList<T, Size>::size() const
{
    return Size;
}


template<class T, unsigned Size>
inline Foam::label Foam::FixedList<T, Size>::max_size() const
{
    return Size;
}


template<class T, unsigned Size>
inline bool Foam::FixedList<T, Size>::empty() const
{
    return false;
}


template<class T, unsigned Size>
template<class HashT>
inline unsigned Foam::FixedList<T, Size>::Hash<HashT>::operator()
(
    const FixedList<T, Size>& lst,
    unsigned seed
) const
{
    if (contiguous<T>())
    {
        // hash directly
        return Hasher(lst.v_, sizeof(lst.v_), seed);
    }
    else
    {
        // hash incrementally
        unsigned val = seed;

        for (register unsigned i=0; i<Size; i++)
        {
            val = HashT()(lst[i], val);
        }

        return val;
    }
}


// ************************************************************************* //
