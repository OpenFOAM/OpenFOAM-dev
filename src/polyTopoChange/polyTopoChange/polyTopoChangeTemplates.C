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

#include "polyTopoChange.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::polyTopoChange::reorder
(
    const labelList& oldToNew,
    DynamicList<T>& lst
)
{
    // Create copy
    DynamicList<T> oldLst(lst);

    forAll(oldToNew, elemI)
    {
        label newElemI = oldToNew[elemI];

        if (newElemI != -1)
        {
            lst[newElemI] = oldLst[elemI];
        }
    }
}


template<class T>
void Foam::polyTopoChange::reorder
(
    const labelList& oldToNew,
    List<DynamicList<T>>& lst
)
{
    // Create copy
    List<DynamicList<T>> oldLst(lst);

    forAll(oldToNew, elemI)
    {
        label newElemI = oldToNew[elemI];

        if (newElemI != -1)
        {
            lst[newElemI].transfer(oldLst[elemI]);
        }
    }
}


template<class T>
void Foam::polyTopoChange::renumberKey
(
    const labelList& oldToNew,
    Map<T>& elems
)
{
    Map<T> newElems(elems.size());

    forAllConstIter(typename Map<T>, elems, iter)
    {
        label newElem = oldToNew[iter.key()];

        if (newElem >= 0)
        {
            newElems.insert(newElem, iter());
        }
    }

    elems.transfer(newElems);
}


// ************************************************************************* //
