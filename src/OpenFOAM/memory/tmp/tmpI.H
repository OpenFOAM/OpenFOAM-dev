/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    isTmp_(true),
    ptr_(tPtr),
    ref_(*tPtr)
{}


template<class T>
inline Foam::tmp<T>::tmp(const T& tRef)
:
    isTmp_(false),
    ptr_(0),
    ref_(tRef)
{}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>& t)
:
    isTmp_(t.isTmp_),
    ptr_(t.ptr_),
    ref_(t.ref_)
{
    if (isTmp_)
    {
        if (ptr_)
        {
            ptr_->operator++();
        }
        else
        {
            FatalErrorIn("Foam::tmp<T>::tmp(const tmp<T>&)")
                << "attempted copy of a deallocated temporary"
                << " of type " << typeid(T).name()
                << abort(FatalError);
        }
    }
}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>& t, bool allowTransfer)
:
    isTmp_(t.isTmp_),
    ptr_(t.ptr_),
    ref_(t.ref_)
{
    if (isTmp_)
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
                FatalErrorIn
                (
                    "Foam::tmp<T>::tmp(const tmp<T>&, bool allowTransfer)"
                )   << "attempted copy of a deallocated temporary"
                    << " of type " << typeid(T).name()
                    << abort(FatalError);
            }
        }
    }
}


template<class T>
inline Foam::tmp<T>::~tmp()
{
    if (isTmp_ && ptr_)
    {
        if (ptr_->okToDelete())
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
    return isTmp_;
}


template<class T>
inline bool Foam::tmp<T>::empty() const
{
    return (isTmp_ && !ptr_);
}


template<class T>
inline bool Foam::tmp<T>::valid() const
{
    return (!isTmp_ || (isTmp_ && ptr_));
}


template<class T>
inline T* Foam::tmp<T>::ptr() const
{
    if (isTmp_)
    {
         if (!ptr_)
         {
             FatalErrorIn("Foam::tmp<T>::ptr() const")
                 << "temporary of type " << typeid(T).name() << " deallocated"
                 << abort(FatalError);
         }

         T* ptr = ptr_;
         ptr_ = 0;

         ptr->resetRefCount();

         return ptr;
    }
    else
    {
        return new T(ref_);
    }
}


template<class T>
inline void Foam::tmp<T>::clear() const
{
    if (isTmp_ && ptr_)  // skip this bit:  && ptr_->okToDelete())
    {
        delete ptr_;
        ptr_ = 0;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline T& Foam::tmp<T>::operator()()
{
    if (isTmp_)
    {
        if (!ptr_)
        {
            FatalErrorIn("T& Foam::tmp<T>::operator()()")
                << "temporary of type " << typeid(T).name() << " deallocated"
                << abort(FatalError);
        }

        return *ptr_;
    }
    else
    {
        // Note: const is cast away!
        // Perhaps there should be two refs, one for const and one for non const
        // and if the ref is actually const then you cannot return it here.
        //
        // Another possibility would be to store a const ref and a flag to say
        // whether the tmp was constructed with a const or a non-const argument.
        //
        // eg, enum refType { POINTER = 0, REF = 1, CONSTREF = 2 };
        return const_cast<T&>(ref_);
    }
}


template<class T>
inline const T& Foam::tmp<T>::operator()() const
{
    if (isTmp_)
    {
        if (!ptr_)
        {
            FatalErrorIn("const T& Foam::tmp<T>::operator()() const")
                << "temporary of type " << typeid(T).name() << " deallocated"
                << abort(FatalError);
        }

        return *ptr_;
    }
    else
    {
        return ref_;
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
    if (isTmp_)
    {
         if (!ptr_)
         {
             FatalErrorIn("Foam::tmp<T>::operator->()")
                 << "temporary of type " << typeid(T).name() << " deallocated"
                 << abort(FatalError);
         }

         return ptr_;
    }
    else
    {
        return &const_cast<T&>(ref_);
    }
}


template<class T>
inline const T* Foam::tmp<T>::operator->() const
{
    return const_cast<tmp<T>&>(*this).operator->();
}


template<class T>
inline void Foam::tmp<T>::operator=(const tmp<T>& t)
{
    if (isTmp_ && ptr_)
    {
        if (ptr_->okToDelete())
        {
            delete ptr_;
            ptr_ = 0;
        }
        else
        {
            ptr_->operator--();
        }
    }

    if (t.isTmp_)
    {
        isTmp_ = true;
        ptr_ = t.ptr_;

        if (ptr_)
        {
            ptr_->operator++();
        }
        else
        {
            FatalErrorIn("Foam::tmp<T>::operator=(const tmp<T>&)")
                << "attempted copy of a deallocated temporary"
                << " of type " << typeid(T).name()
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn("Foam::tmp<T>::operator=(const tmp<T>&)")
            << "attempted to assign to a const reference to constant object"
            << " of type " << typeid(T).name()
            << abort(FatalError);
    }
}


// ************************************************************************* //
