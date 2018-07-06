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

#include "ILList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::ILList<LListBase, T>::ILList(const ILList<LListBase, T>& lst)
:
    UILList<LListBase, T>()
{
    for
    (
        typename UILList<LListBase, T>::const_iterator iter = lst.begin();
        iter != lst.end();
        ++iter
    )
    {
        this->append(iter().clone().ptr());
    }
}


template<class LListBase, class T>
template<class CloneArg>
Foam::ILList<LListBase, T>::ILList
(
    const ILList<LListBase, T>& lst,
    const CloneArg& cloneArg
)
:
    UILList<LListBase, T>()
{
    for
    (
        typename UILList<LListBase, T>::const_iterator iter = lst.begin();
        iter != lst.end();
        ++iter
    )
    {
        this->append(iter().clone(cloneArg).ptr());
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::ILList<LListBase, T>::~ILList()
{
    this->clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class LListBase, class T>
bool Foam::ILList<LListBase, T>::eraseHead()
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
bool Foam::ILList<LListBase, T>::erase(T* p)
{
    T* tPtr;
    if ((tPtr = remove(p)))
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
void Foam::ILList<LListBase, T>::clear()
{
    label oldSize = this->size();
    for (label i=0; i<oldSize; ++i)
    {
        eraseHead();
    }

    LListBase::clear();
}


template<class LListBase, class T>
void Foam::ILList<LListBase, T>::transfer(ILList<LListBase, T>& lst)
{
    clear();
    LListBase::transfer(lst);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class LListBase, class T>
void Foam::ILList<LListBase, T>::operator=(const ILList<LListBase, T>& lst)
{
    this->clear();

    for
    (
        typename UILList<LListBase, T>::const_iterator iter = lst.begin();
        iter != lst.end();
        ++iter
    )
    {
        this->append(iter().clone().ptr());
    }
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "ILListIO.C"


// ************************************************************************* //
