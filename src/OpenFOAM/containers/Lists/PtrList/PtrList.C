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

#include "PtrList.H"
#include "SLPtrList.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::PtrList<T>::PtrList()
:
    ptrs_()
{}


template<class T>
Foam::PtrList<T>::PtrList(const label s)
:
    ptrs_(s, reinterpret_cast<T*>(0))
{}


template<class T>
Foam::PtrList<T>::PtrList(const PtrList<T>& a)
:
    ptrs_(a.size())
{
    forAll(*this, i)
    {
        ptrs_[i] = (a[i]).clone().ptr();
    }
}


template<class T>
template<class CloneArg>
Foam::PtrList<T>::PtrList(const PtrList<T>& a, const CloneArg& cloneArg)
:
    ptrs_(a.size())
{
    forAll(*this, i)
    {
        ptrs_[i] = (a[i]).clone(cloneArg).ptr();
    }
}


template<class T>
Foam::PtrList<T>::PtrList(const Xfer<PtrList<T> >& lst)
{
    transfer(lst());
}


template<class T>
Foam::PtrList<T>::PtrList(PtrList<T>& a, bool reUse)
:
    ptrs_(a.size())
{
    if (reUse)
    {
        forAll(*this, i)
        {
            ptrs_[i] = a.ptrs_[i];
            a.ptrs_[i] = NULL;
        }
        a.setSize(0);
    }
    else
    {
        forAll(*this, i)
        {
            ptrs_[i] = (a[i]).clone().ptr();
        }
    }
}


template<class T>
Foam::PtrList<T>::PtrList(const SLPtrList<T>& sll)
:
    ptrs_(sll.size())
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
            ptrs_[i++] = (iter()).clone().ptr();
        }
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class T>
Foam::PtrList<T>::~PtrList()
{
    forAll(*this, i)
    {
        if (ptrs_[i])
        {
            delete ptrs_[i];
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::PtrList<T>::setSize(const label newSize)
{
    if (newSize < 0)
    {
        FatalErrorIn("PtrList<T>::setSize(const label)")
            << "bad set size " << newSize
            << " for type " << typeid(T).name()
            << abort(FatalError);
    }

    label oldSize = size();

    if (newSize == 0)
    {
        clear();
    }
    else if (newSize < oldSize)
    {
        register label i;
        for (i=newSize; i<oldSize; i++)
        {
            if (ptrs_[i])
            {
                delete ptrs_[i];
            }
        }

        ptrs_.setSize(newSize);
    }
    else // newSize > oldSize
    {
        ptrs_.setSize(newSize);

        register label i;
        for (i=oldSize; i<newSize; i++)
        {
            ptrs_[i] = NULL;
        }
    }
}


template<class T>
void Foam::PtrList<T>::clear()
{
    forAll(*this, i)
    {
        if (ptrs_[i])
        {
            delete ptrs_[i];
        }
    }

    ptrs_.clear();
}


template<class T>
void Foam::PtrList<T>::transfer(PtrList<T>& a)
{
    clear();
    ptrs_.transfer(a.ptrs_);
}


template<class T>
void Foam::PtrList<T>::reorder(const labelUList& oldToNew)
{
    if (oldToNew.size() != size())
    {
        FatalErrorIn("PtrList<T>::reorder(const labelUList&)")
            << "Size of map (" << oldToNew.size()
            << ") not equal to list size (" << size()
            << ") for type " << typeid(T).name()
            << abort(FatalError);
    }

    List<T*> newPtrs_(ptrs_.size(), reinterpret_cast<T*>(0));

    forAll(*this, i)
    {
        label newI = oldToNew[i];

        if (newI < 0 || newI >= size())
        {
            FatalErrorIn("PtrList<T>::reorder(const labelUList&)")
                << "Illegal index " << newI << nl
                << "Valid indices are 0.." << size()-1
                << " for type " << typeid(T).name()
                << abort(FatalError);
        }

        if (newPtrs_[newI])
        {
            FatalErrorIn("PtrList<T>::reorder(const labelUList&)")
                << "reorder map is not unique; element " << newI
                << " already set for type " << typeid(T).name()
                << abort(FatalError);
        }
        newPtrs_[newI] = ptrs_[i];
    }

    forAll(newPtrs_, i)
    {
        if (!newPtrs_[i])
        {
            FatalErrorIn("PtrList<T>::reorder(const labelUList&)")
                << "Element " << i << " not set after reordering with type "
                << typeid(T).name() << nl << abort(FatalError);
        }
    }

    ptrs_.transfer(newPtrs_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
Foam::PtrList<T>& Foam::PtrList<T>::operator=(const PtrList<T>& a)
{
    if (this == &a)
    {
        FatalErrorIn("PtrList<T>::operator=(const PtrList<T>&)")
            << "attempted assignment to self for type " << typeid(T).name()
            << abort(FatalError);
    }

    if (size() == 0)
    {
        setSize(a.size());

        forAll(*this, i)
        {
            ptrs_[i] = (a[i]).clone().ptr();
        }
    }
    else if (a.size() == size())
    {
        forAll(*this, i)
        {
            (*this)[i] = a[i];
        }
    }
    else
    {
        FatalErrorIn("PtrList::operator=(const PtrList<T>&)")
            << "bad size: " << a.size()
            << " for type " << typeid(T).name()
            << abort(FatalError);
    }


    return *this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PtrListIO.C"

// ************************************************************************* //
