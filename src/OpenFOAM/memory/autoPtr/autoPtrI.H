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

#include "error.H"
#include <typeinfo>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::autoPtr<T>::autoPtr(T* p)
:
    ptr_(p)
{}


template<class T>
inline Foam::autoPtr<T>::autoPtr(const autoPtr<T>& ap)
:
    ptr_(ap.ptr_)
{
    ap.ptr_ = nullptr;
}


template<class T>
inline Foam::autoPtr<T>::autoPtr(const autoPtr<T>& ap, const bool reuse)
{
    if (reuse)
    {
        ptr_ = ap.ptr_;
        ap.ptr_ = nullptr;
    }
    else if (ap.valid())
    {
        ptr_ = ap().clone().ptr();
    }
    else
    {
        ptr_ = nullptr;
    }
}


template<class T>
inline Foam::autoPtr<T>::~autoPtr()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::autoPtr<T>::empty() const
{
    return !ptr_;
}


template<class T>
inline bool Foam::autoPtr<T>::valid() const
{
    return ptr_;
}


template<class T>
inline T* Foam::autoPtr<T>::ptr()
{
    T* ptr = ptr_;
    ptr_ = nullptr;
    return ptr;
}


template<class T>
inline void Foam::autoPtr<T>::set(T* p)
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
inline void Foam::autoPtr<T>::reset(T* p)
{
    if (ptr_)
    {
        delete ptr_;
    }

    ptr_ = p;
}


template<class T>
inline void Foam::autoPtr<T>::clear()
{
    reset(nullptr);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline T& Foam::autoPtr<T>::operator()()
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
inline const T& Foam::autoPtr<T>::operator()() const
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
inline T& Foam::autoPtr<T>::operator*()
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
inline const T& Foam::autoPtr<T>::operator*() const
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
inline Foam::autoPtr<T>::operator const T&() const
{
    return operator()();
}


template<class T>
inline T* Foam::autoPtr<T>::operator->()
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
inline const T* Foam::autoPtr<T>::operator->() const
{
    return const_cast<autoPtr<T>&>(*this).operator->();
}


template<class T>
inline void Foam::autoPtr<T>::operator=(T* p)
{
    reset(p);
}


template<class T>
inline void Foam::autoPtr<T>::operator=(const autoPtr<T>& ap)
{
    if (this != &ap)
    {
        reset(const_cast<autoPtr<T>&>(ap).ptr());
    }
}


// ************************************************************************* //
