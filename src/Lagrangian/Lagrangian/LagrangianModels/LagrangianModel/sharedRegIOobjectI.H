/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "sharedRegIOobject.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::sharedRegIOobject<Type>::sharedRegIOobject()
:
    ptr_(nullptr)
{}


template<class Type>
template<class ... Args>
inline Foam::sharedRegIOobject<Type>::sharedRegIOobject
(
    const IOobject& io,
    const Args& ... args
)
:
    ptr_(nullptr)
{
    set(io, args ...);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::sharedRegIOobject<Type>::~sharedRegIOobject()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
inline bool Foam::sharedRegIOobject<Type>::valid() const
{
    return ptr_;
}


template<class Type>
template<class ... Args>
inline void Foam::sharedRegIOobject<Type>::set
(
    const IOobject& io,
    const Args& ... args
)
{
    clear();

    if (!io.db().found(io.name()))
    {
        ptr_ = new Type(io, args ...);

        ptr_->checkIn();
    }
    else
    {
        ptr_ = &io.db().lookupObjectRef<Type>(io.name());

        ptr_->operator++();
    }
}


template<class Type>
inline Type& Foam::sharedRegIOobject<Type>::ref() const
{
    if (!valid())
    {
        FatalErrorInFunction
            << "Attempted access of invalid object"
            << exit(FatalError);
    }

    return *ptr_;
}


template<class Type>
inline void Foam::sharedRegIOobject<Type>::clear()
{
    if (valid())
    {
        if (ptr_->unique())
        {
            ptr_->checkOut();
        }
        else
        {
            ptr_->operator--();
        }

        ptr_ = nullptr;
    }
}


template<class Type>
inline const Type& Foam::sharedRegIOobject<Type>::operator()() const
{
    if (!valid())
    {
        FatalErrorInFunction
            << "Attempted access of invalid object"
            << exit(FatalError);
    }

    return *ptr_;
}


// ************************************************************************* //
