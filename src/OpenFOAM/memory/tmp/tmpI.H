/*--------------------------------------------------------------------r-------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
inline Foam::tmp<T>::tmp(T* tPtr)
:
    type_(TMP),
    ptr_(tPtr)
{}


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
    if (isTmp())
    {
        if (ptr_)
        {
            ptr_->operator++();
        }
        else
        {
            FatalErrorInFunction
                << "Attempted copy of a deallocated temporary"
                << " of type " << typeid(T).name()
                << abort(FatalError);
        }
    }
}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>& t, bool allowTransfer)
:
    type_(t.type_),
    ptr_(t.ptr_)
{
    if (isTmp())
    {
        if (allowTransfer)
        {
            const_cast<tmp<T>&>(t).ptr_ = 0;
        }
        else
        {
            if (ptr_)
            {
                ptr_->operator++();
            }
            else
            {
                FatalErrorInFunction
                    << "Attempted copy of a deallocated temporary"
                    << " of type " << typeid(T).name()
                    << abort(FatalError);
            }
        }
    }
}


template<class T>
inline Foam::tmp<T>::~tmp()
{
    if (isTmp() && ptr_)
    {
        if (ptr_->unique())
        {
            delete ptr_;
            ptr_ = 0;
        }
        else
        {
            ptr_->operator--();
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::tmp<T>::isTmp() const
{
    return type_ == TMP;
}


template<class T>
inline bool Foam::tmp<T>::empty() const
{
    return (isTmp() && !ptr_);
}


template<class T>
inline bool Foam::tmp<T>::valid() const
{
    return (!isTmp() || (isTmp() && ptr_));
}


template<class T>
inline T& Foam::tmp<T>::ref()
{
    if (type_ == TMP)
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << "Temporary of type " << typeid(T).name() << " deallocated"
                << abort(FatalError);
        }

        return *ptr_;
    }
    else
    {
        FatalErrorInFunction
            << "Attempt to acquire non-const reference to const object"
            << abort(FatalError);
        return *ptr_;
    }
}


template<class T>
inline T* Foam::tmp<T>::ptr() const
{
    if (isTmp())
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << "Temporary of type " << typeid(T).name() << " deallocated"
                << abort(FatalError);
        }

        if (!ptr_->unique())
        {
            FatalErrorInFunction
                << "Attempt to acquire pointer to object referred to"
                   " by multiple 'tmp's"
                << abort(FatalError);
        }

        T* ptr = ptr_;
        ptr_ = 0;

        return ptr;
    }
    else
    {
        return new T(*ptr_);
    }
}


template<class T>
inline void Foam::tmp<T>::clear() const
{
    if (isTmp() && ptr_)
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

#ifndef CONST_TMP
template<class T>
inline T& Foam::tmp<T>::operator()()
{
    if (type_ == TMP)
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << "Temporary of type " << typeid(T).name() << " deallocated"
                << abort(FatalError);
        }

        return *ptr_;
    }
    else
    {
        // Const-ness is automatically cast-away which is why this operator is
        // deprecated.  Use ref() where non-const access is required.
        return *ptr_;
    }
}
#endif


template<class T>
inline const T& Foam::tmp<T>::operator()() const
{
    if (type_ == TMP)
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << "Temporary of type " << typeid(T).name() << " deallocated"
                << abort(FatalError);
        }

        return *ptr_;
    }
    else
    {
        // Return const reference
        return *ptr_;
    }
}


template<class T>
inline Foam::tmp<T>::operator const T&() const
{
    return operator()();
}


template<class T>
inline T* Foam::tmp<T>::operator->()
{
    if (isTmp())
    {
         if (!ptr_)
         {
             FatalErrorInFunction
                 << "Temporary of type " << typeid(T).name() << " deallocated"
                 << abort(FatalError);
         }

         return ptr_;
    }
    else
    {
        FatalErrorInFunction << "Const object cast to non-const"
            << abort(FatalError);
        return ptr_;
    }
}


template<class T>
inline const T* Foam::tmp<T>::operator->() const
{
    if (isTmp())
    {
         if (!ptr_)
         {
             FatalErrorInFunction
                 << "Temporary of type " << typeid(T).name() << " deallocated"
                 << abort(FatalError);
         }

         return ptr_;
    }
    else
    {
        return ptr_;
    }
}


template<class T>
inline void Foam::tmp<T>::operator=(T* tPtr)
{
    if (isTmp() && ptr_)
    {
        if (ptr_->unique())
        {
            delete ptr_;
            ptr_ = 0;
        }
        else
        {
            ptr_->operator--();
        }
    }

    type_ = TMP;

    if (!tPtr)
    {
        FatalErrorInFunction
            << "Attempted copy of a deallocated temporary"
            << " of type " << typeid(T).name()
            << abort(FatalError);
    }

    ptr_ = tPtr;
    ptr_->resetRefCount();
}


template<class T>
inline void Foam::tmp<T>::operator=(const tmp<T>& t)
{
    if (isTmp() && ptr_)
    {
        if (ptr_->unique())
        {
            delete ptr_;
            ptr_ = 0;
        }
        else
        {
            ptr_->operator--();
        }
    }

    if (t.isTmp())
    {
        type_ = TMP;

        if (!t.ptr_)
        {
            FatalErrorInFunction
                << "Attempted assignment to a deallocated temporary"
                << " of type " << typeid(T).name()
                << abort(FatalError);
        }

        ptr_ = t.ptr_;

        ptr_->operator++();
    }
    else
    {
        FatalErrorInFunction
            << "Attempted assignment to a const reference to constant object"
            << " of type " << typeid(T).name()
            << abort(FatalError);
    }
}


//- Return the const reference of the non-const reference argument
template<typename T>
inline const T& Const(T& t)
{
    return t;
}


//- Return the const reference of the non-const rvalue reference argument
template<typename T>
inline const T& Const(T&& t)
{
    return t;
}


// ************************************************************************* //
