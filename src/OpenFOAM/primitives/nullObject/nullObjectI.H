/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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

template<class T>
inline const T& Foam::NullObjectRef()
{
    return *reinterpret_cast<const T*>(nullObjectPtr);
}

template<class T>
inline const T* Foam::NullObjectPtr()
{
    return reinterpret_cast<const T*>(nullObjectPtr);
}


template<class T>
inline bool Foam::isNull(const T& t)
{
    return &t == NullObjectPtr<T>();
}

template<class T>
inline bool Foam::notNull(const T& t)
{
    return &t != NullObjectPtr<T>();
}


template<class T>
inline bool Foam::isNull(const T* t)
{
    return t == NullObjectPtr<T>();
}

template<class T>
inline bool Foam::notNull(const T* t)
{
    return t != NullObjectPtr<T>();
}


// ************************************************************************* //
