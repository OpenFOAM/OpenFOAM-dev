/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "patchToPatchTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class SubListA, class SubListB>
void Foam::patchToPatchTools::transferListList
(
    List<SubListA>& a,
    List<SubListB>& b
)
{
    a.resize(b.size());
    forAll(a, i)
    {
        a[i].transfer(b[i]);
    }
}


template<class Type>
void Foam::patchToPatchTools::rDistributeListList
(
    const label size,
    const distributionMap& map,
    List<List<Type>>& data
)
{
    distributionMapBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>(),
        size,
        map.constructMap(),
        false,
        map.subMap(),
        false,
        data,
        ListAppendEqOp<Type>(),
        flipOp(),
        List<Type>()
    );
}


template<class Type>
void Foam::patchToPatchTools::rDistributeListList
(
    const label size,
    const distributionMap& map,
    List<DynamicList<Type>>& data
)
{
    List<List<Type>> tData;
    transferListList(tData, data);
    rDistributeListList(size, map, tData);
    transferListList(data, tData);
}


// ************************************************************************* //
