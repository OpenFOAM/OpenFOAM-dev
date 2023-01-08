/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "ULPtrList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::ULPtrList<LListBase, T>::ULPtrList(const ULPtrList<LListBase, T>& lst)
{
    for (const_iterator iter = lst.begin(); iter != lst.end(); ++iter)
    {
        this->append(&iter());
    }
}


template<class LListBase, class T>
Foam::ULPtrList<LListBase, T>::ULPtrList(ULPtrList<LListBase, T>&& lst)
{
    transfer(lst);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class LListBase, class T>
void Foam::ULPtrList<LListBase, T>::operator=
(
    const ULPtrList<LListBase, T>& lst
)
{
    for (const_iterator iter = lst.begin(); iter != lst.end(); ++iter)
    {
        this->append(&iter());
    }
}


template<class LListBase, class T>
void Foam::ULPtrList<LListBase, T>::operator=(ULPtrList<LListBase, T>&& lst)
{
    transfer(lst);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "ULPtrListIO.C"


// ************************************************************************* //
