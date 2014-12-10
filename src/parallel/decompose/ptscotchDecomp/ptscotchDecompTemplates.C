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

#include "ptscotchDecomp.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Insert at front of list
template<class Type>
void Foam::ptscotchDecomp::prepend
(
    const UList<Type>& extraLst,
    List<Type>& lst
)
{
    label nExtra = extraLst.size();

    // Make space for initial elements
    lst.setSize(lst.size() + nExtra);
    for (label i = lst.size()-1; i >= nExtra; i--)
    {
        lst[i] = lst[i-nExtra];
    }

    // Insert at front
    forAll(extraLst, i)
    {
        lst[i] = extraLst[i];
    }
}


// Insert at back of list
template<class Type>
void Foam::ptscotchDecomp::append
(
    const UList<Type>& extraLst,
    List<Type>& lst
)
{
    label sz = lst.size();

    // Make space for initial elements
    lst.setSize(sz + extraLst.size());

    // Insert at back
    forAll(extraLst, i)
    {
        lst[sz++] = extraLst[i];
    }
}


// ************************************************************************* //
