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

#include "LPtrList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::LPtrList<LListBase, T>::LPtrList(const LPtrList<LListBase, T>& lst)
:
    LList<LListBase, T*>()
{
    for (const_iterator iter = lst.begin(); iter != lst.end(); ++iter)
    {
        this->append(iter().clone().ptr());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::LPtrList<LListBase, T>::~LPtrList()
{
    this->clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class LListBase, class T>
bool Foam::LPtrList<LListBase, T>::eraseHead()
{
    T* tPtr;
    if ((tPtr = this->removeHead()))
    {
        delete tPtr;
        return true;
    }
    else
    {
        return false;
    }
}


template<class LListBase, class T>
void Foam::LPtrList<LListBase, T>::clear()
{
    const label oldSize = this->size();
    for (label i=0; i<oldSize; ++i)
    {
        eraseHead();
    }

    LList<LListBase, T*>::clear();
}


template<class LListBase, class T>
void Foam::LPtrList<LListBase, T>::transfer(LPtrList<LListBase, T>& lst)
{
    clear();
    LList<LListBase, T*>::transfer(lst);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class LListBase, class T>
void Foam::LPtrList<LListBase, T>::operator=(const LPtrList<LListBase, T>& lst)
{
    clear();

    for (const_iterator iter = lst.begin(); iter != lst.end(); ++iter)
    {
        this->append(iter().clone().ptr());
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "LPtrListIO.C"


// ************************************************************************* //
