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

#include "List.H"
#include "ListLoopM.H"

#include "FixedList.H"
#include "PtrList.H"
#include "SLList.H"
#include "IndirectList.H"
#include "UIndirectList.H"
#include "BiIndirectList.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

// Construct with length specified
template<class T>
Foam::List<T>::List(const label s)
:
    UList<T>(NULL, s)
{
    if (this->size_ < 0)
    {
        FatalErrorIn("List<T>::List(const label size)")
            << "bad size " << this->size_
            << abort(FatalError);
    }

    if (this->size_)
    {
        this->v_ = new T[this->size_];
    }
}


// Construct with length and single value specified
template<class T>
Foam::List<T>::List(const label s, const T& a)
:
    UList<T>(NULL, s)
{
    if (this->size_ < 0)
    {
        FatalErrorIn("List<T>::List(const label size, const T&)")
            << "bad size " << this->size_
            << abort(FatalError);
    }

    if (this->size_)
    {
        this->v_ = new T[this->size_];

        List_ACCESS(T, (*this), vp);
        List_FOR_ALL((*this), i)
            List_ELEM((*this), vp, i) = a;
        List_END_FOR_ALL
    }
}


// Construct as copy
template<class T>
Foam::List<T>::List(const List<T>& a)
:
    UList<T>(NULL, a.size_)
{
    if (this->size_)
    {
        this->v_ = new T[this->size_];

#       ifdef USEMEMCPY
        if (contiguous<T>())
        {
            memcpy(this->v_, a.v_, this->byteSize());
        }
        else
#       endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            List_FOR_ALL((*this), i)
                List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
            List_END_FOR_ALL
        }
    }
}


// Construct by transferring the parameter contents
template<class T>
Foam::List<T>::List(const Xfer<List<T> >& lst)
{
    transfer(lst());
}


// Construct as copy or re-use as specified.
template<class T>
Foam::List<T>::List(List<T>& a, bool reUse)
:
    UList<T>(NULL, a.size_)
{
    if (reUse)
    {
        this->v_ = a.v_;
        a.v_ = 0;
        a.size_ = 0;
    }
    else if (this->size_)
    {
        this->v_ = new T[this->size_];

#       ifdef USEMEMCPY
        if (contiguous<T>())
        {
            memcpy(this->v_, a.v_, this->byteSize());
        }
        else
#       endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            List_FOR_ALL((*this), i)
                List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
            List_END_FOR_ALL
        }
    }
}


// Construct as subset
template<class T>
Foam::List<T>::List(const UList<T>& a, const labelUList& map)
:
    UList<T>(NULL, map.size())
{
    if (this->size_)
    {
        // Note:cannot use List_ELEM since third argument has to be index.

        this->v_ = new T[this->size_];

        forAll(*this, i)
        {
            this->v_[i] = a[map[i]];
        }
    }
}


// Construct given start and end iterators.
template<class T>
template<class InputIterator>
Foam::List<T>::List(InputIterator first, InputIterator last)
{
    label s = 0;
    for
    (
        InputIterator iter = first;
        iter != last;
        ++iter
    )
    {
        s++;
    }

    setSize(s);

    s = 0;

    for
    (
        InputIterator iter = first;
        iter != last;
        ++iter
    )
    {
        this->operator[](s++) = *iter;
    }
}


// Construct as copy of FixedList<T, Size>
template<class T>
template<unsigned Size>
Foam::List<T>::List(const FixedList<T, Size>& lst)
:
    UList<T>(NULL, Size)
{
    if (this->size_)
    {
        this->v_ = new T[this->size_];

        forAll(*this, i)
        {
            this->operator[](i) = lst[i];
        }
    }
}


// Construct as copy of PtrList<T>
template<class T>
Foam::List<T>::List(const PtrList<T>& lst)
:
    UList<T>(NULL, lst.size())
{
    if (this->size_)
    {
        this->v_ = new T[this->size_];

        forAll(*this, i)
        {
            this->operator[](i) = lst[i];
        }
    }
}


// Construct as copy of SLList<T>
template<class T>
Foam::List<T>::List(const SLList<T>& lst)
:
    UList<T>(NULL, lst.size())
{
    if (this->size_)
    {
        this->v_ = new T[this->size_];

        label i = 0;
        for
        (
            typename SLList<T>::const_iterator iter = lst.begin();
            iter != lst.end();
            ++iter
        )
        {
            this->operator[](i++) = iter();
        }
    }
}


// Construct as copy of UIndirectList<T>
template<class T>
Foam::List<T>::List(const UIndirectList<T>& lst)
:
    UList<T>(NULL, lst.size())
{
    if (this->size_)
    {
        this->v_ = new T[this->size_];

        forAll(*this, i)
        {
            this->operator[](i) = lst[i];
        }
    }
}


// Construct as copy of BiIndirectList<T>
template<class T>
Foam::List<T>::List(const BiIndirectList<T>& lst)
:
    UList<T>(NULL, lst.size())
{
    if (this->size_)
    {
        this->v_ = new T[this->size_];

        forAll(*this, i)
        {
            this->operator[](i) = lst[i];
        }
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

// Destroy list elements
template<class T>
Foam::List<T>::~List()
{
    if (this->v_) delete[] this->v_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::List<T>::setSize(const label newSize)
{
    if (newSize < 0)
    {
        FatalErrorIn("List<T>::setSize(const label)")
            << "bad set size " << newSize
            << abort(FatalError);
    }

    if (newSize != this->size_)
    {
        if (newSize > 0)
        {
            T* nv = new T[label(newSize)];

            if (this->size_)
            {
                register label i = min(this->size_, newSize);

#               ifdef USEMEMCPY
                if (contiguous<T>())
                {
                    memcpy(nv, this->v_, i*sizeof(T));
                }
                else
#               endif
                {
                    register T* vv = &this->v_[i];
                    register T* av = &nv[i];
                    while (i--) *--av = *--vv;
                }
            }
            if (this->v_) delete[] this->v_;

            this->size_ = newSize;
            this->v_ = nv;
        }
        else
        {
            clear();
        }
    }
}


template<class T>
void Foam::List<T>::setSize(const label newSize, const T& a)
{
    label oldSize = label(this->size_);
    this->setSize(newSize);

    if (newSize > oldSize)
    {
        register label i = newSize - oldSize;
        register T* vv = &this->v_[newSize];
        while (i--) *--vv = a;
    }
}


template<class T>
void Foam::List<T>::clear()
{
    if (this->v_) delete[] this->v_;
    this->size_ = 0;
    this->v_ = 0;
}


// Transfer the contents of the argument List into this List
// and annul the argument list
template<class T>
void Foam::List<T>::transfer(List<T>& a)
{
    if (this->v_) delete[] this->v_;
    this->size_ = a.size_;
    this->v_ = a.v_;

    a.size_ = 0;
    a.v_ = 0;
}


// Transfer the contents of the argument DynamicList into this List
// and annul the argument list
template<class T>
template<unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
void Foam::List<T>::transfer(DynamicList<T, SizeInc, SizeMult, SizeDiv>& a)
{
    // shrink the allocated space to the number of elements used
    a.shrink();
    transfer(static_cast<List<T>&>(a));
    a.clearStorage();
}


// Transfer the contents of the argument SortableList into this List
// and annul the argument list
template<class T>
void Foam::List<T>::transfer(SortableList<T>& a)
{
    // shrink away the sort indices
    a.shrink();
    transfer(static_cast<List<T>&>(a));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Assignment to UList operator. Takes linear time.
template<class T>
void Foam::List<T>::operator=(const UList<T>& a)
{
    if (a.size_ != this->size_)
    {
        if (this->v_) delete[] this->v_;
        this->v_ = 0;
        this->size_ = a.size_;
        if (this->size_) this->v_ = new T[this->size_];
    }

    if (this->size_)
    {
#       ifdef USEMEMCPY
        if (contiguous<T>())
        {
            memcpy(this->v_, a.v_, this->byteSize());
        }
        else
#       endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            List_FOR_ALL((*this), i)
                List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
            List_END_FOR_ALL
        }
    }
}


// Assignment operator. Takes linear time.
template<class T>
void Foam::List<T>::operator=(const List<T>& a)
{
    if (this == &a)
    {
        FatalErrorIn("List<T>::operator=(const List<T>&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    operator=(static_cast<const UList<T>&>(a));
}


// Assignment operator. Takes linear time.
template<class T>
void Foam::List<T>::operator=(const SLList<T>& lst)
{
    if (lst.size() != this->size_)
    {
        if (this->v_) delete[] this->v_;
        this->v_ = 0;
        this->size_ = lst.size();
        if (this->size_) this->v_ = new T[this->size_];
    }

    if (this->size_)
    {
        label i = 0;
        for
        (
            typename SLList<T>::const_iterator iter = lst.begin();
            iter != lst.end();
            ++iter
        )
        {
            this->operator[](i++) = iter();
        }
    }
}


// Assignment operator. Takes linear time.
template<class T>
void Foam::List<T>::operator=(const UIndirectList<T>& lst)
{
    if (lst.size() != this->size_)
    {
        if (this->v_) delete[] this->v_;
        this->v_ = 0;
        this->size_ = lst.size();
        if (this->size_) this->v_ = new T[this->size_];
    }

    forAll(*this, i)
    {
        this->operator[](i) = lst[i];
    }
}


// Assignment operator. Takes linear time.
template<class T>
void Foam::List<T>::operator=(const BiIndirectList<T>& lst)
{
    if (lst.size() != this->size_)
    {
        if (this->v_) delete[] this->v_;
        this->v_ = 0;
        this->size_ = lst.size();
        if (this->size_) this->v_ = new T[this->size_];
    }

    forAll(*this, i)
    {
        this->operator[](i) = lst[i];
    }
}

// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "ListIO.C"

// ************************************************************************* //
