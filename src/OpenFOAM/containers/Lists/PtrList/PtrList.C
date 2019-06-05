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

#include "PtrList.H"
#include "SLPtrList.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::PtrList<T>::PtrList()
:
    UPtrList<T>()
{}


template<class T>
Foam::PtrList<T>::PtrList(const label s)
:
    UPtrList<T>(s)
{}


template<class T>
Foam::PtrList<T>::PtrList(const PtrList<T>& a)
:
    UPtrList<T>(a.size())
{
    forAll(*this, i)
    {
        this->ptrs_[i] = (a[i]).clone().ptr();
    }
}


template<class T>
template<class CloneArg>
Foam::PtrList<T>::PtrList(const PtrList<T>& a, const CloneArg& cloneArg)
:
    UPtrList<T>(a.size())
{
    forAll(*this, i)
    {
        this->ptrs_[i] = (a[i]).clone(cloneArg).ptr();
    }
}


template<class T>
Foam::PtrList<T>::PtrList(PtrList<T>&& lst)
{
    transfer(lst);
}


template<class T>
Foam::PtrList<T>::PtrList(PtrList<T>& a, bool reuse)
:
    UPtrList<T>(a, reuse)
{
    if (!reuse)
    {
        forAll(*this, i)
        {
            this->ptrs_[i] = (a[i]).clone().ptr();
        }
    }
}


template<class T>
Foam::PtrList<T>::PtrList(const SLPtrList<T>& sll)
:
    UPtrList<T>(sll.size())
{
    if (sll.size())
    {
        label i = 0;
        for
        (
            typename SLPtrList<T>::const_iterator iter = sll.begin();
            iter != sll.end();
            ++iter
        )
        {
            this->ptrs_[i++] = (iter()).clone().ptr();
        }
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class T>
Foam::PtrList<T>::~PtrList()
{
    forAll(*this, i)
    {
        if (this->ptrs_[i])
        {
            delete this->ptrs_[i];
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::PtrList<T>::setSize(const label newSize)
{
    if (newSize < 0)
    {
        FatalErrorInFunction
            << "bad set size " << newSize
            << " for type " << typeid(T).name()
            << abort(FatalError);
    }

    label oldSize = this->size();

    if (newSize == 0)
    {
        clear();
    }
    else if (newSize < oldSize)
    {
        label i;
        for (i=newSize; i<oldSize; i++)
        {
            if (this->ptrs_[i])
            {
                delete this->ptrs_[i];
            }
        }

        this->ptrs_.setSize(newSize);
    }
    else // newSize > oldSize
    {
        this->ptrs_.setSize(newSize);

        label i;
        for (i=oldSize; i<newSize; i++)
        {
            this->ptrs_[i] = nullptr;
        }
    }
}


template<class T>
void Foam::PtrList<T>::clear()
{
    forAll(*this, i)
    {
        if (this->ptrs_[i])
        {
            delete this->ptrs_[i];
        }
    }

    this->ptrs_.clear();
}


template<class T>
void Foam::PtrList<T>::transfer(PtrList<T>& a)
{
    clear();
    this->ptrs_.transfer(a.ptrs_);
}


template<class T>
void Foam::PtrList<T>::reorder(const labelUList& oldToNew)
{
    if (oldToNew.size() != this->size())
    {
        FatalErrorInFunction
            << "Size of map (" << oldToNew.size()
            << ") not equal to list size (" << this->size()
            << ") for type " << typeid(T).name()
            << abort(FatalError);
    }

    List<T*> newPtrs_(this->ptrs_.size(), reinterpret_cast<T*>(0));

    forAll(*this, i)
    {
        label newI = oldToNew[i];

        if (newI < 0 || newI >= this->size())
        {
            FatalErrorInFunction
                << "Illegal index " << newI << nl
                << "Valid indices are 0.." << this->size()-1
                << " for type " << typeid(T).name()
                << abort(FatalError);
        }

        if (newPtrs_[newI])
        {
            FatalErrorInFunction
                << "reorder map is not unique; element " << newI
                << " already set for type " << typeid(T).name()
                << abort(FatalError);
        }
        newPtrs_[newI] = this->ptrs_[i];
    }

    forAll(newPtrs_, i)
    {
        if (!newPtrs_[i])
        {
            FatalErrorInFunction
                << "Element " << i << " not set after reordering with type "
                << typeid(T).name() << nl << abort(FatalError);
        }
    }

    this->ptrs_.transfer(newPtrs_);
}


template<class T>
void Foam::PtrList<T>::shuffle(const labelUList& newToOld)
{
    List<T*> newPtrs_(newToOld.size(), reinterpret_cast<T*>(0));

    forAll(newToOld, newI)
    {
        label oldI = newToOld[newI];

        if (oldI >= 0 && oldI < this->size())
        {
            newPtrs_[newI] = this->ptrs_[oldI];
            this->ptrs_[oldI] = nullptr;
        }
    }

    // Delete all remaining pointers
    clear();

    // Take over new pointers
    this->ptrs_.transfer(newPtrs_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::PtrList<T>::operator=(const PtrList<T>& a)
{
    if (this == &a)
    {
        FatalErrorInFunction
            << "attempted assignment to self for type " << typeid(T).name()
            << abort(FatalError);
    }

    if (this->size() == 0)
    {
        setSize(a.size());

        forAll(*this, i)
        {
            this->ptrs_[i] = (a[i]).clone().ptr();
        }
    }
    else if (a.size() == this->size())
    {
        forAll(*this, i)
        {
            (*this)[i] = a[i];
        }
    }
    else
    {
        FatalErrorInFunction
            << "bad size: " << a.size()
            << " for type " << typeid(T).name()
            << abort(FatalError);
    }
}


template<class T>
void Foam::PtrList<T>::operator=(PtrList<T>&& a)
{
    if (this == &a)
    {
        FatalErrorInFunction
            << "attempted assignment to self for type " << typeid(T).name()
            << abort(FatalError);
    }

    transfer(a);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PtrListIO.C"

// ************************************************************************* //
