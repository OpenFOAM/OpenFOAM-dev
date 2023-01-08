/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "UPtrListDictionary.H"
#include "PtrListDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::UPtrListDictionary<T>::UPtrListDictionary(const label size)
:
    DictionaryBase<UPtrList<T>, T>(2*size)
{
    UPtrList<T>::setSize(size);
}


template<class T>
Foam::UPtrListDictionary<T>::UPtrListDictionary(const UPtrListDictionary& dict)
:
    DictionaryBase<UPtrList<T>, T>(dict)
{}


template<class T>
Foam::UPtrListDictionary<T>::UPtrListDictionary
(
    const PtrListDictionary<T>& dict
)
:
    DictionaryBase<UPtrList<T>, T>(0)
{
    DictionaryBase<UPtrList<T>, T>::hashedTs_ = dict.hashedTs_;
    UPtrList<T>::operator=(dict);
}


template<class T>
Foam::UPtrListDictionary<T>::UPtrListDictionary(UPtrListDictionary&& dict)
:
    DictionaryBase<UPtrList<T>, T>(move(dict))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::UPtrListDictionary<T>::set(const label i) const
{
    return UPtrList<T>::set(i);
}


template<class T>
inline T* Foam::UPtrListDictionary<T>::set
(
    const label i,
    const word& key,
    T* ptr
)
{
    if (!DictionaryBase<UPtrList<T>, T>::hashedTs_.insert(key, ptr))
    {
        FatalErrorInFunction
            << "Cannot insert with key '" << key << "' into hash-table"
            << abort(FatalError);
    }
    return UPtrList<T>::set(i, ptr);
}


// ************************************************************************* //
