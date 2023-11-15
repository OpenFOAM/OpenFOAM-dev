/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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


template<class Type, class LabelList, class ScalarList>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatchTools::interpolate
(
    const List<LabelList>& localOtherFaces,
    const List<ScalarList>& weights,
    const autoPtr<distributionMap>& otherMapPtr,
    const Field<Type>& otherFld
)
{
    // Distribute the other field if necessary
    tmp<Field<Type>> tLocalOtherFld;
    if (otherMapPtr.valid())
    {
        tLocalOtherFld = tmp<Field<Type>>(new Field<Type>(otherFld));
        otherMapPtr->distribute(tLocalOtherFld.ref());
    }
    const Field<Type>& localOtherFld =
        otherMapPtr.valid() ? tLocalOtherFld() : otherFld;

    // Allocate the result
    tmp<Field<Type>> tFld
    (
        new Field<Type>(localOtherFaces.size(), pTraits<Type>::nan)
    );
    Field<Type>& fld = tFld.ref();

    // Compute the result as a weighted sum
    forAll(localOtherFaces, facei)
    {
        scalar sumW = 0;
        Type sumWF = Zero;

        forAll(localOtherFaces[facei], i)
        {
            const scalar w = weights[facei][i];
            sumW += w;
            sumWF += w*localOtherFld[localOtherFaces[facei][i]];
        }

        if (localOtherFaces[facei].size())
        {
            fld[facei] = sumWF/sumW;
        }
    }

    return tFld;
}


template<class Type, class LabelList, class ScalarList>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatchTools::interpolate
(
    const List<LabelList>& localOtherFaces,
    const List<ScalarList>& weights,
    const autoPtr<distributionMap>& otherMapPtr,
    const Field<Type>& otherFld,
    const Field<Type>& leftOverFld
)
{
    // Distribute the other field if necessary
    tmp<Field<Type>> tLocalOtherFld;
    if (otherMapPtr.valid())
    {
        tLocalOtherFld = tmp<Field<Type>>(new Field<Type>(otherFld));
        otherMapPtr->distribute(tLocalOtherFld.ref());
    }
    const Field<Type>& localOtherFld =
        otherMapPtr.valid() ? tLocalOtherFld() : otherFld;

    // Allocate the result
    tmp<Field<Type>> tFld
    (
        new Field<Type>(localOtherFaces.size(), pTraits<Type>::nan)
    );
    Field<Type>& fld = tFld.ref();

    // Compute the result as a weighted sum
    forAll(localOtherFaces, facei)
    {
        scalar sumW = 0;
        Type sumWF = Zero;

        forAll(localOtherFaces[facei], i)
        {
            const scalar w = weights[facei][i];
            sumW += w;
            sumWF += w*localOtherFld[localOtherFaces[facei][i]];
        }

        fld[facei] = sumWF + (1 - sumW)*leftOverFld[facei];
    }

    return tFld;
}


// ************************************************************************* //
