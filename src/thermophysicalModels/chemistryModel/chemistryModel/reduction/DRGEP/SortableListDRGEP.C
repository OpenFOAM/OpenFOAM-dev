/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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

#include "OSspecific.H"
#include <algorithm>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Type>
Foam::SortableListDRGEP<Type>::SortableListDRGEP(const List<Type>& values)
:
    List<Type>(values),
    indices_(values.size())
{
    sort();
}


template <class Type>
Foam::SortableListDRGEP<Type>::SortableListDRGEP(const label size)
:
    List<Type>(size),
    indices_(size)
{}


template <class Type>
Foam::SortableListDRGEP<Type>::SortableListDRGEP
(
    const label size,
    const Type& val
)
:
    List<Type>(size, val),
    indices_(size)
{}


template <class Type>
Foam::SortableListDRGEP<Type>::SortableListDRGEP
(
    const SortableListDRGEP<Type>& lst
)
:
    List<Type>(lst),
    indices_(lst.indices())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class Type>
void Foam::SortableListDRGEP<Type>::setSize(const label newSize)
{
    List<Type>::setSize(newSize);
    indices_.setSize(newSize);
}


template <class Type>
void Foam::SortableListDRGEP<Type>::sort()
{
    forAll(indices_, i)
    {
        indices_[i] = i;
    }

    Foam::sort(indices_, less(*this));

    List<Type> tmpValues(this->size());

    forAll(indices_, i)
    {
        tmpValues[i] = this->operator[](indices_[i]);
    }

    List<Type>::transfer(tmpValues);
}


template <class Type>
void Foam::SortableListDRGEP<Type>::partialSort(int M)
{
    forAll(indices_, i)
    {
        indices_[i] = i;
    }

    std::partial_sort
    (
        indices_.begin(),
        indices_.begin()+M,
        indices_.end(),
        less(*this)
    );
}


template <class Type>
void Foam::SortableListDRGEP<Type>::stableSort()
{
    forAll(indices_, i)
    {
        indices_[i] = i;
    }

    Foam::stableSort(indices_, less(*this));

    List<Type> tmpValues(this->size());

    forAll(indices_, i)
    {
        tmpValues[i] = this->operator[](indices_[i]);
    }

    List<Type>::transfer(tmpValues);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template <class Type>
void
Foam::SortableListDRGEP<Type>::operator=(const SortableListDRGEP<Type>& rhs)
{
    List<Type>::operator=(rhs);
    indices_ = rhs.indices();
}


// ************************************************************************* //
