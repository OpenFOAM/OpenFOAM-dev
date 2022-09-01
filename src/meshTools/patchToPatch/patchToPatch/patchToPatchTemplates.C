/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "patchToPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Type Foam::patchToPatch::NaN()
{
    Type result;

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        setComponent(result, cmpt) = Foam::NaN;
    }

    return result;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatch::interpolate
(
    const List<DynamicList<label>>& localOtherFaces,
    const List<DynamicList<scalar>>& weights,
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
    tmp<Field<Type>> tFld(new Field<Type>(localOtherFaces.size(), NaN<Type>()));
    Field<Type>& fld = tFld.ref();

    // Compute the result as a weighted sum
    forAll(localOtherFaces, tgtFacei)
    {
        scalar sumW = 0;
        Type sumWF = Zero;

        forAll(localOtherFaces[tgtFacei], i)
        {
            const scalar w = weights[tgtFacei][i];
            sumW += w;
            sumWF += w*localOtherFld[localOtherFaces[tgtFacei][i]];
        }

        if (localOtherFaces[tgtFacei].size())
        {
            fld[tgtFacei] = sumWF/sumW;
        }
    }

    return tFld;
}

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatch::interpolate
(
    const List<DynamicList<label>>& localOtherFaces,
    const List<DynamicList<scalar>>& weights,
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
    tmp<Field<Type>> tFld(new Field<Type>(localOtherFaces.size(), NaN<Type>()));
    Field<Type>& fld = tFld.ref();

    // Compute the result as a weighted sum
    forAll(localOtherFaces, tgtFacei)
    {
        scalar sumW = 0;
        Type sumWF = Zero;

        forAll(localOtherFaces[tgtFacei], i)
        {
            const scalar w = weights[tgtFacei][i];
            sumW += w;
            sumWF += w*localOtherFld[localOtherFaces[tgtFacei][i]];
        }

        fld[tgtFacei] = sumWF + (1 - sumW)*leftOverFld[tgtFacei];
    }

    return tFld;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchToPatch::srcToTgt(const Field<Type>& srcFld) const
{
    return
        interpolate
        (
            tgtLocalSrcFaces_,
            tgtWeights(),
            srcMapPtr_,
            srcFld
        );
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatch::srcToTgt
(
    const Field<Type>& srcFld,
    const Field<Type>& leftOverTgtFld
) const
{
    return
        interpolate
        (
            tgtLocalSrcFaces_,
            tgtWeights(),
            srcMapPtr_,
            srcFld,
            leftOverTgtFld
        );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchToPatch::tgtToSrc(const Field<Type>& tgtFld) const
{
    return
        interpolate
        (
            srcLocalTgtFaces_,
            srcWeights(),
            tgtMapPtr_,
            tgtFld
        );
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatch::tgtToSrc
(
    const Field<Type>& tgtFld,
    const Field<Type>& leftOverSrcFld
) const
{
    return
        interpolate
        (
            srcLocalTgtFaces_,
            srcWeights(),
            tgtMapPtr_,
            tgtFld,
            leftOverSrcFld
        );
}
// ************************************************************************* //
