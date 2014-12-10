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

Description
    Patch to patch interpolation functions

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Interpolate point field
template<class FromPatch, class ToPatch>
template<class Type>
tmp<Field<Type> >
PatchToPatchInterpolation<FromPatch, ToPatch>::pointInterpolate
(
    const Field<Type>& pf
) const
{
    if (pf.size() != fromPatch_.nPoints())
    {
        FatalErrorIn
        (
            "PatchToPatchInterpolation::pointInterpolate"
            "(const Field<Type> pf)"
        )   << "given field does not correspond to patch. Patch size: "
            << fromPatch_.nPoints() << " field size: " << pf.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            toPatch_.nPoints(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    const List<typename FromPatch::FaceType>& fromPatchLocalFaces =
        fromPatch_.localFaces();

    const FieldField<Field, scalar>& weights = pointWeights();

    const labelList& addr = pointAddr();

    forAll(result, pointI)
    {
        const scalarField& curWeights = weights[pointI];

        if (addr[pointI] > -1)
        {
            const labelList& hitFacePoints =
                fromPatchLocalFaces[addr[pointI]];

            forAll(curWeights, wI)
            {
                result[pointI] += curWeights[wI]*pf[hitFacePoints[wI]];
            }
        }
    }

    return tresult;
}


template<class FromPatch, class ToPatch>
template<class Type>
tmp<Field<Type> >
PatchToPatchInterpolation<FromPatch, ToPatch>::pointInterpolate
(
    const tmp<Field<Type> >& tpf
) const
{
    tmp<Field<Type> > tint = pointInterpolate<Type>(tpf());
    tpf.clear();
    return tint;
}


//- Interpolate face field
template<class FromPatch, class ToPatch>
template<class Type>
tmp<Field<Type> >
PatchToPatchInterpolation<FromPatch, ToPatch>::faceInterpolate
(
    const Field<Type>& ff
) const
{
    if (ff.size() != fromPatch_.size())
    {
        FatalErrorIn
        (
            "PatchToPatchInterpolation::faceInterpolate"
            "(const Field<Type> ff)"
        )   << "given field does not correspond to patch. Patch size: "
            << fromPatch_.size() << " field size: " << ff.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            toPatch_.size(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    const labelListList& fromPatchFaceFaces = fromPatch_.faceFaces();

    const FieldField<Field, scalar>& weights = faceWeights();

    const labelList& addr = faceAddr();

    forAll(result, faceI)
    {
        const scalarField& curWeights = weights[faceI];

        if (addr[faceI] > -1)
        {
            const labelList& hitFaceFaces =
                fromPatchFaceFaces[addr[faceI]];

            // first add the hit face
            result[faceI] += ff[addr[faceI]]*curWeights[0];

            for (label wI = 1; wI < curWeights.size(); wI++)
            {
                result[faceI] += ff[hitFaceFaces[wI - 1]]*curWeights[wI];
            }
        }
    }

    return tresult;
}


template<class FromPatch, class ToPatch>
template<class Type>
tmp<Field<Type> >
PatchToPatchInterpolation<FromPatch, ToPatch>::faceInterpolate
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = faceInterpolate(tff());
    tff.clear();
    return tint;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
