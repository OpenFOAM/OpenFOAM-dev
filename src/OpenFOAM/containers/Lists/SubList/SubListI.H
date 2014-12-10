/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::SubList<T>::SubList
(
    const UList<T>& list,
    const label subSize
)
:
    UList<T>(list.v_, subSize)
{
#   ifdef FULLDEBUG
    list.checkSize(subSize);
#   endif
}


template<class T>
inline Foam::SubList<T>::SubList
(
    const UList<T>& list,
    const label subSize,
    const label startIndex
)
:
    UList<T>(&(list.v_[startIndex]), subSize)
{
#   ifdef FULLDEBUG

    // Artificially allow the start of a zero-sized subList to be
    // one past the end of the original list.
    if (subSize)
    {
        list.checkStart(startIndex);
        list.checkSize(startIndex + subSize);
    }
    else
    {
        // Start index needs to fall between 0 and size.  One position
        // behind the last element is allowed
        list.checkSize(startIndex);
    }
#   endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const Foam::SubList<T>& Foam::SubList<T>::null()
{
    return *reinterpret_cast< SubList<T>* >(0);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline Foam::SubList<T>::operator const Foam::List<T>&() const
{
    return *reinterpret_cast< const List<T>* >(this);
}


template<class T>
inline void Foam::SubList<T>::operator=(const T& t)
{
    UList<T>::operator=(t);
}


// ************************************************************************* //
