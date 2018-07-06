/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "ListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::SortableList<T>::SortableList()
{}


template<class T>
Foam::SortableList<T>::SortableList(const UList<T>& values)
:
    List<T>(values)
{
    sort();
}


template<class T>
Foam::SortableList<T>::SortableList(const Xfer<List<T>>& values)
:
    List<T>(values)
{
    sort();
}


template<class T>
Foam::SortableList<T>::SortableList(const label size)
:
    List<T>(size)
{}


template<class T>
Foam::SortableList<T>::SortableList(const label size, const T& val)
:
    List<T>(size, val)
{}


template<class T>
Foam::SortableList<T>::SortableList(const SortableList<T>& lst)
:
    List<T>(lst),
    indices_(lst.indices())
{}


template<class T>
Foam::SortableList<T>::SortableList(std::initializer_list<T> values)
:
    List<T>(values)
{
    sort();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class T>
void Foam::SortableList<T>::clear()
{
    List<T>::clear();
    indices_.clear();
}


template<class T>
Foam::List<T>& Foam::SortableList<T>::shrink()
{
    indices_.clear();
    return static_cast<List<T>&>(*this);
}


template<class T>
void Foam::SortableList<T>::sort()
{
    sortedOrder(*this, indices_);

    List<T> lst(this->size());
    forAll(indices_, i)
    {
        lst[i] = this->operator[](indices_[i]);
    }

    List<T>::transfer(lst);
}


template<class T>
void Foam::SortableList<T>::reverseSort()
{
    sortedOrder(*this, indices_, typename UList<T>::greater(*this));

    List<T> lst(this->size());
    forAll(indices_, i)
    {
        lst[i] = this->operator[](indices_[i]);
    }

    List<T>::transfer(lst);
}


template<class T>
Foam::Xfer<Foam::List<T>> Foam::SortableList<T>::xfer()
{
    return xferMoveTo<List<T>>(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::SortableList<T>::operator=(const T& t)
{
    UList<T>::operator=(t);
}


template<class T>
inline void Foam::SortableList<T>::operator=(const UList<T>& lst)
{
    List<T>::operator=(lst);
    indices_.clear();
}


template<class T>
inline void Foam::SortableList<T>::operator=(const SortableList<T>& lst)
{
    List<T>::operator=(lst);
    indices_ = lst.indices();
}


template<class T>
inline void Foam::SortableList<T>::operator=(std::initializer_list<T> lst)
{
    List<T>::operator=(lst);
    sort();
}


// ************************************************************************* //
