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

#include "nullObject.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template<class T>
inline const Foam::Xfer<T>& Foam::Xfer<T>::null()
{
    return NullObjectRef<Xfer<T>>();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::Xfer<T>::Xfer(T* p)
:
    ptr_(p ? p : new T)
{}


template<class T>
inline Foam::Xfer<T>::Xfer(T& t, bool allowTransfer)
:
    ptr_(new T)
{
    if (allowTransfer)
    {
        ptr_->transfer(t);
    }
    else
    {
        ptr_->operator=(t);
    }
}


template<class T>
inline Foam::Xfer<T>::Xfer(const T& t)
:
    ptr_(new T)
{
    ptr_->operator=(t);
}


template<class T>
inline Foam::Xfer<T>::Xfer(const Xfer<T>& t)
:
    ptr_(new T)
{
    ptr_->transfer(*(t.ptr_));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T>
inline Foam::Xfer<T>::~Xfer()
{
    delete ptr_;
    ptr_ = 0;
}


// * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::Xfer<T>::operator=(T& t)
{
    ptr_->transfer(t);
}


template<class T>
inline void Foam::Xfer<T>::operator=(const Xfer<T>& t)
{
    // silently ignore attempted copy to self
    if (this != &t)
    {
        ptr_->transfer(*(t.ptr_));
    }
}


template<class T>
inline T& Foam::Xfer<T>::operator()() const
{
    return *ptr_;
}


template<class T>
inline T* Foam::Xfer<T>::operator->() const
{
    return ptr_;
}


// * * * * * * * * * * * * *  Helper Functions * * * * * * * * * * * * * * * //


template<class T>
inline Foam::Xfer<T> Foam::xferCopy(const T& t)
{
    return Foam::Xfer<T>(t);
}


template<class T>
inline Foam::Xfer<T> Foam::xferMove(T& t)
{
    return Foam::Xfer<T>(t, true);
}


template<class T>
inline Foam::Xfer<T> Foam::xferTmp(Foam::tmp<T>& tt)
{
    return Foam::Xfer<T>(tt(), tt.isTmp());
}


template<class To, class From>
inline Foam::Xfer<To> Foam::xferCopyTo(const From& t)
{
    Foam::Xfer<To> xf;
    xf() = t;
    return xf;
}


template<class To, class From>
inline Foam::Xfer<To> Foam::xferMoveTo(From& t)
{
    Foam::Xfer<To> xf;
    xf().transfer(t);
    return xf;
}


// ************************************************************************* //
