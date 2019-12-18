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

#include "UPtrList.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::UPtrList<T>::UPtrList()
:
    ptrs_()
{}


template<class T>
Foam::UPtrList<T>::UPtrList(const label s)
:
    ptrs_(s, reinterpret_cast<T*>(0))
{}


template<class T>
Foam::UPtrList<T>::UPtrList(UPtrList<T>& a, bool reuse)
:
    ptrs_(a.ptrs_, reuse)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::UPtrList<T>::setSize(const label newSize)
{
    label oldSize = size();

    if (newSize <= 0)
    {
        clear();
    }
    else if (newSize < oldSize)
    {
        ptrs_.setSize(newSize);
    }
    else if (newSize > oldSize)
    {
        ptrs_.setSize(newSize);

        // set new elements to nullptr
        for (label i=oldSize; i<newSize; i++)
        {
            ptrs_[i] = nullptr;
        }
    }
}


template<class T>
void Foam::UPtrList<T>::clear()
{
    ptrs_.clear();
}


template<class T>
void Foam::UPtrList<T>::transfer(UPtrList<T>& a)
{
    ptrs_.transfer(a.ptrs_);
}


template<class T>
void Foam::UPtrList<T>::reorder(const labelUList& oldToNew)
{
    if (oldToNew.size() != size())
    {
        FatalErrorInFunction
            << "Size of map (" << oldToNew.size()
            << ") not equal to list size (" << size()
            << ")." << abort(FatalError);
    }

    List<T*> newPtrs_(ptrs_.size(), reinterpret_cast<T*>(0));

    forAll(*this, i)
    {
        label newI = oldToNew[i];

        if (newI < 0 || newI >= size())
        {
            FatalErrorInFunction
                << "Illegal index " << newI << nl
                << "Valid indices are 0.." << size()-1
                << abort(FatalError);
        }

        if (newPtrs_[newI])
        {
            FatalErrorInFunction
                << "reorder map is not unique; element " << newI
                << " already set." << abort(FatalError);
        }
        newPtrs_[newI] = ptrs_[i];
    }

    forAll(newPtrs_, i)
    {
        if (!newPtrs_[i])
        {
            FatalErrorInFunction
                << "Element " << i << " not set after reordering." << nl
                << abort(FatalError);
        }
    }

    ptrs_.transfer(newPtrs_);
}


template<class T>
void Foam::UPtrList<T>::shuffle(const labelUList& newToOld)
{
    List<T*> newPtrs_(newToOld.size(), reinterpret_cast<T*>(0));

    forAll(newToOld, newI)
    {
        label oldI = newToOld[newI];

        if (oldI >= 0 && oldI < this->size())
        {
            newPtrs_[newI] = this->ptrs_[oldI];
        }
    }

    this->ptrs_.transfer(newPtrs_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UPtrListIO.C"

// ************************************************************************* //
