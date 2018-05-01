/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "Pstream.H"
#include "PstreamBuffers.H"
#include "PstreamCombineReduceOps.H"
#include "globalIndexAndTransform.H"
#include "transformField.H"
#include "flipOp.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
void Foam::mapDistribute::applyDummyTransforms(List<T>& field) const
{
    forAll(transformElements_, trafoI)
    {
        const labelList& elems = transformElements_[trafoI];

        label n = transformStart_[trafoI];

        forAll(elems, i)
        {
            field[n++] = field[elems[i]];
        }
    }
}


template<class T>
void Foam::mapDistribute::applyDummyInverseTransforms(List<T>& field) const
{
    forAll(transformElements_, trafoI)
    {
        const labelList& elems = transformElements_[trafoI];
        label n = transformStart_[trafoI];

        forAll(elems, i)
        {
            field[elems[i]] = field[n++];
        }
    }
}


template<class T, class TransformOp>   //, class CombineOp>
void Foam::mapDistribute::applyTransforms
(
    const globalIndexAndTransform& globalTransforms,
    List<T>& field,
    const TransformOp& top
) const
{
    const List<vectorTensorTransform>& totalTransform =
        globalTransforms.transformPermutations();

    forAll(totalTransform, trafoI)
    {
        const vectorTensorTransform& vt = totalTransform[trafoI];
        const labelList& elems = transformElements_[trafoI];
        label n = transformStart_[trafoI];

        // Could be optimised to avoid memory allocations
        List<T> transformFld(UIndirectList<T>(field, elems));
        top(vt, true, transformFld);

        forAll(transformFld, i)
        {
            // cop(field[n++], transformFld[i]);
            field[n++] = transformFld[i];
        }
    }
}


template<class T, class TransformOp>   //, class CombineOp>
void Foam::mapDistribute::applyInverseTransforms
(
    const globalIndexAndTransform& globalTransforms,
    List<T>& field,
    const TransformOp& top
) const
{
    const List<vectorTensorTransform>& totalTransform =
        globalTransforms.transformPermutations();

    forAll(totalTransform, trafoI)
    {
        const vectorTensorTransform& vt = totalTransform[trafoI];
        const labelList& elems = transformElements_[trafoI];
        label n = transformStart_[trafoI];

        // Could be optimised to avoid memory allocations
        List<T> transformFld(SubList<T>(field, elems.size(), n));
        top(vt, false, transformFld);

        forAll(transformFld, i)
        {
            // cop(field[elems[i]], transformFld[i]);
            field[elems[i]] = transformFld[i];
        }
    }
}


template<class T, class negateOp>
void Foam::mapDistribute::distribute
(
    List<T>& fld,
    const negateOp& negOp,
    const bool dummyTransform,
    const int tag
) const
{
    mapDistributeBase::distribute(fld, negOp, tag);

    //- Fill in transformed slots with copies
    if (dummyTransform)
    {
        applyDummyTransforms(fld);
    }
}


template<class T>
void Foam::mapDistribute::distribute
(
    List<T>& fld,
    const bool dummyTransform,
    const int tag
) const
{
    distribute(fld, flipOp(), dummyTransform, tag);
}


template<class T>
void Foam::mapDistribute::distribute
(
    DynamicList<T>& fld,
    const bool dummyTransform,
    const int tag
) const
{
    fld.shrink();

    List<T>& fldList = static_cast<List<T>& >(fld);

    distribute(fldList, dummyTransform, tag);

    fld.setCapacity(fldList.size());
}


template<class T>
void Foam::mapDistribute::reverseDistribute
(
    const label constructSize,
    List<T>& fld,
    const bool dummyTransform,
    const int tag
) const
{
    if (dummyTransform)
    {
        applyDummyInverseTransforms(fld);
    }

    mapDistributeBase::reverseDistribute(constructSize, fld, tag);
}


template<class T>
void Foam::mapDistribute::reverseDistribute
(
    const label constructSize,
    const T& nullValue,
    List<T>& fld,
    const bool dummyTransform,
    const int tag
) const
{
    if (dummyTransform)
    {
        applyDummyInverseTransforms(fld);
    }

    mapDistributeBase::reverseDistribute(constructSize, nullValue, fld, tag);
}


template<class T, class TransformOp>
void Foam::mapDistribute::distribute
(
    const globalIndexAndTransform& git,
    List<T>& fld,
    const TransformOp& top,
    const int tag
) const
{
    // Distribute. Leave out dummy transforms since we're doing them ourselves
    distribute(fld, false, tag);
    // Do transforms
    applyTransforms(git, fld, top);
}


template<class T, class TransformOp>
void Foam::mapDistribute::reverseDistribute
(
    const globalIndexAndTransform& git,
    const label constructSize,
    List<T>& fld,
    const TransformOp& top,
    const int tag
) const
{
    // Fill slots with reverse-transformed data. Note that it also copies
    // back into the non-remote part of fld even though these values are not
    // used.
    applyInverseTransforms(git, fld, top);

    // And send back (the remote slots). Disable dummy transformations.
    reverseDistribute(constructSize, fld, false, tag);
}


template<class T, class TransformOp>
void Foam::mapDistribute::reverseDistribute
(
    const globalIndexAndTransform& git,
    const label constructSize,
    const T& nullValue,
    List<T>& fld,
    const TransformOp& top,
    const int tag
) const
{
    // Fill slots with reverse-transformed data Note that it also copies
    // back into the non-remote part of fld even though these values are not
    // used.
    applyInverseTransforms(git, fld, top);   //, eqOp<T>());

    // And send back (the remote slots) Disable dummy transformations.
    reverseDistribute(constructSize, nullValue, fld, false, tag);
}


// ************************************************************************* //
