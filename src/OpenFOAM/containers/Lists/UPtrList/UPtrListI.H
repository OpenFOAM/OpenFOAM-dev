/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "error.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::label Foam::UPtrList<T>::size() const
{
    return ptrs_.size();
}


template<class T>
inline bool Foam::UPtrList<T>::empty() const
{
    return ptrs_.empty();
}


template<class T>
inline T& Foam::UPtrList<T>::first()
{
    return this->operator[](0);
}


template<class T>
inline const T& Foam::UPtrList<T>::first() const
{
    return this->operator[](0);
}


template<class T>
inline T& Foam::UPtrList<T>::last()
{
    return this->operator[](this->size()-1);
}


template<class T>
inline const T& Foam::UPtrList<T>::last() const
{
    return this->operator[](this->size()-1);
}


template<class T>
inline void Foam::UPtrList<T>::resize(const label newSize)
{
    this->setSize(newSize);
}


template<class T>
inline bool Foam::UPtrList<T>::set(const label i) const
{
    return ptrs_[i] != NULL;
}


template<class T>
inline T* Foam::UPtrList<T>::set(const label i, T* ptr)
{
    T* old = ptrs_[i];
    ptrs_[i] = ptr;
    return old;
}


template<class T>
inline Foam::Xfer<Foam::UPtrList<T> > Foam::UPtrList<T>::xfer()
{
    return xferMove(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline const T& Foam::UPtrList<T>::operator[](const label i) const
{
    if (!ptrs_[i])
    {
        FatalErrorIn("UPtrList::operator[] const")
            << "hanging pointer at index " << i
            << " (size " << size()
            << "), cannot dereference"
            << abort(FatalError);
    }

    return *(ptrs_[i]);
}


template<class T>
inline T& Foam::UPtrList<T>::operator[](const label i)
{
    if (!ptrs_[i])
    {
        FatalErrorIn("UPtrList::operator[]")
            << "hanging pointer at index " << i
            << " (size " << size()
            << "), cannot dereference"
            << abort(FatalError);
    }

    return *(ptrs_[i]);
}


template<class T>
inline const T* Foam::UPtrList<T>::operator()(const label i) const
{
    return ptrs_[i];
}


// * * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * //

template<class T>
inline Foam::UPtrList<T>::iterator::iterator(T** ptr)
:
    ptr_(ptr)
{}

template<class T>
inline bool Foam::UPtrList<T>::iterator::operator==(const iterator& iter) const
{
    return ptr_ == iter.ptr_;
}

template<class T>
inline bool Foam::UPtrList<T>::iterator::operator!=(const iterator& iter) const
{
    return ptr_ != iter.ptr_;
}

template<class T>
inline T& Foam::UPtrList<T>::iterator::operator*()
{
    return **ptr_;
}

template<class T>
inline T& Foam::UPtrList<T>::iterator::operator()()
{
    return operator*();
}

template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::iterator::operator++()
{
    ++ptr_;
    return *this;
}

template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::iterator::operator++(int)
{
    iterator tmp = *this;
    ++ptr_;
    return tmp;
}

template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::iterator::operator--()
{
    --ptr_;
    return *this;
}

template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::iterator::operator--(int)
{
    iterator tmp = *this;
    --ptr_;
    return tmp;
}

template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::iterator::operator+=(label n)
{
    ptr_ += n;
    return *this;
}

template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::operator+(const typename UPtrList<T>::iterator& iter, label n)
{
    typename UPtrList<T>::iterator tmp = iter;
    return tmp += n;
}

template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::operator+(label n, const typename UPtrList<T>::iterator& iter)
{
    typename UPtrList<T>::iterator tmp = iter;
    return tmp += n;
}

template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::iterator::operator-=(label n)
{
    ptr_ -= n;
    return *this;
}

template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::operator-(const typename UPtrList<T>::iterator& iter, label n)
{
    typename UPtrList<T>::iterator tmp = iter;
    return tmp -= n;
}

template<class T>
inline Foam::label Foam::operator-
(
    const typename UPtrList<T>::iterator& iter1,
    const typename UPtrList<T>::iterator& iter2
)
{
    return (iter1.ptr_ - iter2.ptr_)/sizeof(T*);
}

template<class T>
inline T& Foam::UPtrList<T>::iterator::operator[](label n)
{
    return *(*this + n);
}

template<class T>
inline bool Foam::UPtrList<T>::iterator::operator<(const iterator& iter) const
{
    return ptr_ < iter.ptr_;
}

template<class T>
inline bool Foam::UPtrList<T>::iterator::operator>(const iterator& iter) const
{
    return ptr_ > iter.ptr_;
}

template<class T>
inline bool Foam::UPtrList<T>::iterator::operator<=(const iterator& iter) const
{
    return ptr_ <= iter.ptr_;
}

template<class T>
inline bool Foam::UPtrList<T>::iterator::operator>=(const iterator& iter) const
{
    return ptr_ >= iter.ptr_;
}

template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::begin()
{
    return ptrs_.begin();
}

template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::end()
{
    return ptrs_.end();
}


// ************************************************************************* //
