/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
#include <typeinfo>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::UautoPtr<T>::UautoPtr(T* p)
:
    ptr_(p)
{}


template<class T>
inline Foam::UautoPtr<T>::UautoPtr(const UautoPtr<T>& ap)
:
    ptr_(ap.ptr_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::UautoPtr<T>::empty() const
{
    return !ptr_;
}


template<class T>
inline bool Foam::UautoPtr<T>::valid() const
{
    return ptr_;
}


template<class T>
inline T* Foam::UautoPtr<T>::ptr()
{
    return ptr_;
}


template<class T>
inline void Foam::UautoPtr<T>::set(T* p)
{
    if (ptr_)
    {
        FatalErrorInFunction
            << "object of type " << typeid(T).name()
            << " already allocated"
            << abort(FatalError);
    }

    ptr_ = p;
}


template<class T>
inline void Foam::UautoPtr<T>::reset(T* p)
{
    ptr_ = p;
}


template<class T>
inline void Foam::UautoPtr<T>::clear()
{
    reset(nullptr);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline T& Foam::UautoPtr<T>::operator()()
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << "object of type " << typeid(T).name()
            << " is not allocated"
            << abort(FatalError);
    }

    return *ptr_;
}


template<class T>
inline const T& Foam::UautoPtr<T>::operator()() const
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << "object of type " << typeid(T).name()
            << " is not allocated"
            << abort(FatalError);
    }

    return *ptr_;
}


template<class T>
inline T& Foam::UautoPtr<T>::operator*()
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << "object of type " << typeid(T).name()
            << " is not allocated"
            << abort(FatalError);
    }

    return *ptr_;
}


template<class T>
inline const T& Foam::UautoPtr<T>::operator*() const
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << "object of type " << typeid(T).name()
            << " is not allocated"
            << abort(FatalError);
    }

    return *ptr_;
}


template<class T>
inline Foam::UautoPtr<T>::operator const T&() const
{
    return operator()();
}


template<class T>
inline T* Foam::UautoPtr<T>::operator->()
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << "object of type " << typeid(T).name()
            << " is not allocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline const T* Foam::UautoPtr<T>::operator->() const
{
    return const_cast<UautoPtr<T>&>(*this).operator->();
}


template<class T>
inline void Foam::UautoPtr<T>::operator=(T* p)
{
    reset(p);
}


template<class T>
inline void Foam::UautoPtr<T>::operator=(const UautoPtr<T>& ap)
{
    if (this != &ap)
    {
        reset(const_cast<UautoPtr<T>&>(ap).ptr());
    }
}


// ************************************************************************* //
