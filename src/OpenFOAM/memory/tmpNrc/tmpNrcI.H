/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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
inline Foam::tmpNrc<T>::tmpNrc(T* tPtr)
:
    type_(TMP),
    ptr_(tPtr)
{}


template<class T>
inline Foam::tmpNrc<T>::tmpNrc(const T& tRef)
:
    type_(CONST_REF),
    ptr_(const_cast<T*>(&tRef))
{}


template<class T>
inline Foam::tmpNrc<T>::tmpNrc(const tmpNrc<T>& t)
:
    type_(t.type_),
    ptr_(t.ptr_)
{
    if (isTmp())
    {
        if (ptr_)
        {
            t.type_ = CONST_REF;
        }
        else
        {
            FatalErrorInFunction
                << "Attempted copy of a deallocated " << typeName()
                << abort(FatalError);
        }
    }
}


template<class T>
inline Foam::tmpNrc<T>::tmpNrc(const tmpNrc<T>&& t)
:
    type_(t.type_),
    ptr_(t.ptr_)
{
    if (isTmp())
    {
        t.ptr_ = 0;
    }
}


template<class T>
inline Foam::tmpNrc<T>::tmpNrc(const tmpNrc<T>& t, bool allowTransfer)
:
    type_(t.type_),
    ptr_(t.ptr_)
{
    if (isTmp())
    {
        if (ptr_)
        {
            if (allowTransfer)
            {
                t.ptr_ = 0;
            }
            else
            {
                t.type_ = CONST_REF;
            }
        }
        else
        {
            FatalErrorInFunction
                << "Attempted copy of a deallocated " << typeName()
                << abort(FatalError);
        }
    }
}


template<class T>
inline Foam::tmpNrc<T>::~tmpNrc()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::tmpNrc<T>::isTmp() const
{
    return type_ == TMP;
}


template<class T>
inline bool Foam::tmpNrc<T>::empty() const
{
    return (isTmp() && !ptr_);
}


template<class T>
inline bool Foam::tmpNrc<T>::valid() const
{
    return (!isTmp() || (isTmp() && ptr_));
}


template<class T>
inline Foam::word Foam::tmpNrc<T>::typeName() const
{
    return "tmpNrc<" + word(typeid(T).name()) + '>';
}


template<class T>
inline T& Foam::tmpNrc<T>::ref()
{
    if (isTmp())
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << typeName() << " deallocated"
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Attempt to acquire non-const reference to const object"
            << " from a " << typeName()
            << abort(FatalError);
    }

    return *ptr_;
}


template<class T>
inline T* Foam::tmpNrc<T>::ptr() const
{
    if (isTmp())
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << typeName() << " deallocated"
                << abort(FatalError);
        }

        T* ptr = ptr_;
        ptr_ = 0;

        return ptr;
    }
    else
    {
        return ptr_->clone().ptr();
    }
}


template<class T>
inline void Foam::tmpNrc<T>::clear() const
{
    if (isTmp() && ptr_)
    {
        delete ptr_;
        ptr_ = 0;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline const T& Foam::tmpNrc<T>::operator()() const
{
    if (isTmp())
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << typeName() << " deallocated"
                << abort(FatalError);
        }
    }

    // Return const reference
    return *ptr_;
}


template<class T>
inline Foam::tmpNrc<T>::operator const T&() const
{
    return operator()();
}


template<class T>
inline T* Foam::tmpNrc<T>::operator->()
{
    if (isTmp())
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << typeName() << " deallocated"
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Attempt to cast const object to non-const for a " << typeName()
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline const T* Foam::tmpNrc<T>::operator->() const
{
    if (isTmp() && !ptr_)
    {
        FatalErrorInFunction
            << typeName() << " deallocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline void Foam::tmpNrc<T>::operator=(T* tPtr)
{
    clear();

    if (!tPtr)
    {
        FatalErrorInFunction
            << "Attempted copy of a deallocated " << typeName()
            << abort(FatalError);
    }

    type_ = TMP;
    ptr_ = tPtr;
}


template<class T>
inline void Foam::tmpNrc<T>::operator=(const tmpNrc<T>& t)
{
    clear();

    if (t.isTmp())
    {
        type_ = TMP;

        // if (!t.ptr_)
        // {
        //     FatalErrorInFunction
        //         << "Attempted assignment to a deallocated " << typeName()
        //         << abort(FatalError);
        // }

        ptr_ = t.ptr_;
        t.ptr_ = 0;
    }
    else
    {
        FatalErrorInFunction
            << "Attempted assignment to a const reference to an object"
            << " of type " << typeid(T).name()
            << abort(FatalError);
    }
}


template<class T>
inline void Foam::tmpNrc<T>::operator=(const tmpNrc<T>&& t)
{
    clear();

    type_ = t.type_;
    ptr_ = t.ptr_;

    if (isTmp())
    {
        t.ptr_ = 0;
    }
}


// ************************************************************************* //
