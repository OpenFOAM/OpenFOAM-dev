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

#include "autoPtr.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline void Foam::PtrList<T>::resize(const label newSize)
{
    this->setSize(newSize);
}


template<class T>
inline void Foam::PtrList<T>::append(T* ptr)
{
    label sz = this->size();
    this->setSize(sz+1);
    this->ptrs_[sz] = ptr;
}


template<class T>
inline void Foam::PtrList<T>::append(const autoPtr<T>& aptr)
{
    return append(const_cast<autoPtr<T>&>(aptr).ptr());
}


template<class T>
inline void Foam::PtrList<T>::append
(
    const tmp<T>& t
)
{
    return append(const_cast<tmp<T>&>(t).ptr());
}


template<class T>
inline bool Foam::PtrList<T>::set(const label i) const
{
    return this->ptrs_[i] != nullptr;
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set(const label i, T* ptr)
{
    autoPtr<T> old(this->ptrs_[i]);
    this->ptrs_[i] = ptr;
    return old;
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set
(
    const label i,
    const autoPtr<T>& aptr
)
{
    return set(i, const_cast<autoPtr<T>&>(aptr).ptr());
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set
(
    const label i,
    const tmp<T>& t
)
{
    return set(i, const_cast<tmp<T>&>(t).ptr());
}


// ************************************************************************* //
