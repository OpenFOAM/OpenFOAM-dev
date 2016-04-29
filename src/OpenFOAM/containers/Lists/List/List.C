/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::List<T>::List(const label s)
:
    UList<T>(NULL, s)
{
    if (this->size_ < 0)
    {
        FatalErrorInFunction
            << "bad size " << this->size_
            << abort(FatalError);
    }

    if (this->size_)
    {
        this->v_ = new T[this->size_];
    }
}


template<class T>
Foam::List<T>::List(const label s, const T& a)
:
    UList<T>(NULL, s)
{
    if (this->size_ < 0)
    {
        FatalErrorInFunction
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


template<class T>
Foam::List<T>::List(const label s, const zero)
:
    UList<T>(NULL, s)
{
    if (this->size_ < 0)
    {
        FatalErrorInFunction
            << "bad size " << this->size_
            << abort(FatalError);
    }

    if (this->size_)
    {
        this->v_ = new T[this->size_];

        List_ACCESS(T, (*this), vp);
        List_FOR_ALL((*this), i)
            List_ELEM((*this), vp, i) = Zero;
        List_END_FOR_ALL
    }
}


template<class T>
Foam::List<T>::List(const List<T>& a)
:
    UList<T>(NULL, a.size_)
{
    if (this->size_)
    {
        this->v_ = new T[this->size_];

        #ifdef USEMEMCPY
        if (contiguous<T>())
        {
            memcpy(this->v_, a.v_, this->byteSize());
        }
        else
        #endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            List_FOR_ALL((*this), i)
                List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
            List_END_FOR_ALL
        }
    }
}


template<class T>
Foam::List<T>::List(const Xfer<List<T>>& lst)
{
    transfer(lst());
}


template<class T>
Foam::List<T>::List(List<T>& a, bool reuse)
:
    UList<T>(NULL, a.size_)
{
    if (reuse)
    {
        this->v_ = a.v_;
        a.v_ = 0;
        a.size_ = 0;
    }
    else if (this->size_)
    {
        this->v_ = new T[this->size_];

        #ifdef USEMEMCPY
        if (contiguous<T>())
        {
            memcpy(this->v_, a.v_, this->byteSize());
        }
        else
        #endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            List_FOR_ALL((*this), i)
                List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
            List_END_FOR_ALL
        }
    }
}


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
        FatalErrorInFunction
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
                label i = min(this->size_, newSize);

                #ifdef USEMEMCPY
                if (contiguous<T>())
                {
                    memcpy(nv, this->v_, i*sizeof(T));
                }
                else
                #endif
                {
                    T* vv = &this->v_[i];
                    T* av = &nv[i];
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
        label i = newSize - oldSize;
        T* vv = &this->v_[newSize];
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


template<class T>
void Foam::List<T>::transfer(List<T>& a)
{
    if (this->v_) delete[] this->v_;
    this->size_ = a.size_;
    this->v_ = a.v_;

    a.size_ = 0;
    a.v_ = 0;
}


template<class T>
template<unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
void Foam::List<T>::transfer(DynamicList<T, SizeInc, SizeMult, SizeDiv>& a)
{
    // Shrink the allocated space to the number of elements used
    a.shrink();
    transfer(static_cast<List<T>&>(a));
    a.clearStorage();
}


template<class T>
void Foam::List<T>::transfer(SortableList<T>& a)
{
    // Shrink away the sort indices
    a.shrink();
    transfer(static_cast<List<T>&>(a));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

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
        #ifdef USEMEMCPY
        if (contiguous<T>())
        {
            memcpy(this->v_, a.v_, this->byteSize());
        }
        else
        #endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            List_FOR_ALL((*this), i)
                List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
            List_END_FOR_ALL
        }
    }
}


template<class T>
void Foam::List<T>::operator=(const List<T>& a)
{
    if (this == &a)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    operator=(static_cast<const UList<T>&>(a));
}


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
