/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Operators  * * * * * * * * * * * //

template<class T>
inline void Foam::tmp<T>::operator++()
{
    ptr_->operator++();

    if (ptr_->count() > 1)
    {
        FatalErrorInFunction
            << "Attempt to create more than 2 tmp's referring to"
               " the same object of type " << typeName()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::tmp<T>::tmp(T* tPtr, bool nonReusable)
:
    type_(nonReusable ? NON_RESUSABLE_TMP : REUSABLE_TMP),
    ptr_(tPtr)
{
    if (tPtr && !tPtr->unique())
    {
        FatalErrorInFunction
            << "Attempted construction of a " << typeName()
            << " from non-unique pointer"
            << abort(FatalError);
    }
}


template<class T>
inline Foam::tmp<T>::tmp(const T& tRef)
:
    type_(CONST_REF),
    ptr_(const_cast<T*>(&tRef))
{}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>& t)
:
    type_(t.type_),
    ptr_(t.ptr_)
{
    if (isAnyTmp())
    {
        if (ptr_)
        {
            operator++();
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
inline Foam::tmp<T>::tmp(const tmp<T>&& t)
:
    type_(t.type_),
    ptr_(t.ptr_)
{
    if (isAnyTmp())
    {
        t.ptr_ = 0;
    }
}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>& t, bool allowTransfer)
:
    type_(t.type_),
    ptr_(t.ptr_)
{
    if (isAnyTmp())
    {
        if (ptr_)
        {
            if (allowTransfer && type_ == REUSABLE_TMP)
            {
                t.ptr_ = 0;
            }
            else
            {
                operator++();
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
inline Foam::tmp<T>::~tmp()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::tmp<T>::isAnyTmp() const
{
    return type_ == REUSABLE_TMP || type_ == NON_RESUSABLE_TMP;
}


template<class T>
inline bool Foam::tmp<T>::isTmp() const
{
    return type_ == REUSABLE_TMP;
}


template<class T>
inline bool Foam::tmp<T>::empty() const
{
    return (isAnyTmp() && !ptr_);
}


template<class T>
inline bool Foam::tmp<T>::valid() const
{
    return (!isAnyTmp() || (isAnyTmp() && ptr_));
}


template<class T>
inline Foam::word Foam::tmp<T>::typeName() const
{
    return "tmp<" + word(typeid(T).name()) + '>';
}


template<class T>
inline T& Foam::tmp<T>::ref() const
{
    if (isAnyTmp())
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
inline T* Foam::tmp<T>::ptr() const
{
    if (isTmp())
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << typeName() << " deallocated"
                << abort(FatalError);
        }

        if (!ptr_->unique())
        {
            FatalErrorInFunction
                << "Attempt to acquire pointer to object referred to"
                << " by multiple temporaries of type " << typeName()
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
inline void Foam::tmp<T>::clear() const
{
    if (isAnyTmp() && ptr_)
    {
        if (ptr_->unique())
        {
            delete ptr_;
            ptr_ = 0;
        }
        else
        {
            ptr_->operator--();
            ptr_ = 0;
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

#ifdef NON_CONST_TMP
template<class T>
inline T& Foam::tmp<T>::operator()()
{
    if (isAnyTmp())
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << typeName() << " deallocated"
                << abort(FatalError);
        }
    }

    // Const-ness is automatically cast-away which is why this operator is
    // deprecated.  Use ref() where non-const access is required.
    return *ptr_;
}
#endif


template<class T>
inline const T& Foam::tmp<T>::operator()() const
{
    if (isAnyTmp())
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
inline Foam::tmp<T>::operator const T&() const
{
    return operator()();
}


template<class T>
inline T* Foam::tmp<T>::operator->()
{
    if (isAnyTmp())
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
inline const T* Foam::tmp<T>::operator->() const
{
    if (isAnyTmp() && !ptr_)
    {
        FatalErrorInFunction
            << typeName() << " deallocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline void Foam::tmp<T>::operator=(T* tPtr)
{
    clear();

    if (!tPtr)
    {
        FatalErrorInFunction
            << "Attempted copy of a deallocated " << typeName()
            << abort(FatalError);
    }

    if (tPtr && !tPtr->unique())
    {
        FatalErrorInFunction
            << "Attempted assignment of a " << typeName()
            << " to non-unique pointer"
            << abort(FatalError);
    }

    type_ = REUSABLE_TMP;
    ptr_ = tPtr;
}


template<class T>
inline void Foam::tmp<T>::operator=(const tmp<T>& t)
{
    clear();

    if (t.isAnyTmp())
    {
        type_ = t.type_;

        if (!t.ptr_)
        {
            FatalErrorInFunction
                << "Attempted assignment to a deallocated " << typeName()
                << abort(FatalError);
        }

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


// ************************************************************************* //
