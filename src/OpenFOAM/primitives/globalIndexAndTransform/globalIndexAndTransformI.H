/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "polyMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::globalIndexAndTransform::less::operator()
(
    const labelPair& a,
    const labelPair& b
) const
{
    label procA = globalIndexAndTransform::processor(a);
    label procB = globalIndexAndTransform::processor(b);

    if (procA < procB)
    {
        return true;
    }
    else if (procA > procB)
    {
        return false;
    }
    else
    {
        // Equal proc.
        label indexA = globalIndexAndTransform::index(a);
        label indexB = globalIndexAndTransform::index(b);

        if (indexA < indexB)
        {
            return true;
        }
        else if (indexA > indexB)
        {
            return false;
        }
        else
        {
            // Equal index
            label transformA = globalIndexAndTransform::transformIndex(a);
            label transformB = globalIndexAndTransform::transformIndex(b);

            return transformA < transformB;
        }
    }
}


Foam::label Foam::globalIndexAndTransform::encodeTransformIndex
(
    const List<label>& permutationIndices
) const
{
    if (permutationIndices.size() != transforms_.size())
    {
        FatalErrorIn
        (
            "Foam::label encodeTransformIndex"
            "("
                "const List<label>& permutationIndices,"
            ") const"
        )
            << "permutationIndices " << permutationIndices
            << "are of a different size to the number of independent transforms"
            << abort(FatalError);
    }

    label transformIndex = 0;

    label w = 1;

    forAll(transforms_, b)
    {
        if (mag(permutationIndices[b]) > 1)
        {
            FatalErrorIn
            (
                "Foam::label encodeTransformIndex"
                "("
                "const List<label>& permutationIndices,"
                ") const"
            )
                << "permutationIndices " << permutationIndices
                << "are illegal, they must all be only -1, 0 or +1"
                << abort(FatalError);
        }

        transformIndex += (permutationIndices[b] + 1)*w;

        w *= 3;
    }

    return transformIndex;
}


Foam::label Foam::globalIndexAndTransform::encodeTransformIndex
(
    const FixedList<Foam::label, 3>& permutation
) const
{
    if (nIndependentTransforms() == 0)
    {
        return 0;
    }
    if (nIndependentTransforms() == 1)
    {
        return permutation[0]+1;
    }
    else if (nIndependentTransforms() == 2)
    {
        return (permutation[1]+1)*3 + (permutation[0]+1);
    }
    else
    {
        return
            (permutation[2]+1)*9
          + (permutation[1]+1)*3
          + (permutation[0]+1);
    }
}


Foam::FixedList<Foam::label, 3>
Foam::globalIndexAndTransform::decodeTransformIndex
(
    const label transformIndex
) const
{
    FixedList<label, 3> permutation(label(0));

    label t = transformIndex;
    if (nIndependentTransforms() > 0)
    {
        permutation[0] = (t%3)-1;
        if (nIndependentTransforms() > 1)
        {
            t /= 3;
            permutation[1] = (t%3)-1;
            if (nIndependentTransforms() > 2)
            {
                t /= 3;
                permutation[2] = (t%3)-1;
            }
        }
    }

#   ifdef FULLDEBUG
    t /= 3;
    if (t != 0)
    {
        FatalErrorIn
        (
            "globalIndexAndTransform::decodeTransformIndex(const label)"
        )   << "transformIndex : " << transformIndex
            << " has more than 3 fields."
            << abort(FatalError);
    }
#   endif

    return permutation;
}


Foam::label Foam::globalIndexAndTransform::addToTransformIndex
(
    const label transformIndex,
    const label patchI,
    const bool isSendingSide,
    const scalar tol
) const
{
    const Pair<label>& transSign = patchTransformSign_[patchI];

    label matchTransI = transSign.first();

    // Hardcoded for max 3 transforms only!

    if (matchTransI > -1 && matchTransI < 3)
    {
        FixedList<label, 3> permutation = decodeTransformIndex(transformIndex);


        // Add patch transform
        // ~~~~~~~~~~~~~~~~~~~

        label sign = transSign.second();
        if (!isSendingSide)
        {
            sign = -sign;
        }


        // If this transform been found already by a patch?
        if (permutation[matchTransI] != 0)
        {
            if (sign == 0)
            {
                // sent from patch without a transformation. Do nothing.
                FatalErrorIn("globalIndexAndTransform::addToTransformIndex(..)")
                    << "patch:" << mesh_.boundaryMesh()[patchI].name()
                    << " transform:" << matchTransI << " sign:" << sign
                    << "  current transforms:" << permutation
                    << exit(FatalError);
            }
            else if (sign == permutation[matchTransI])
            {
                // This is usually illegal. The only exception is for points
                // on the axis of a 180 degree cyclic wedge when the
                // transformation is going to be (-1 0 0 0 -1 0 0 0 +1)
                // (or a different permutation but always two times -1 and
                // once +1)
                bool antiCyclic = false;

                const vectorTensorTransform& vt = transforms_[matchTransI];
                if (mag(vt.t()) < SMALL && vt.hasR())
                {
                    const tensor& R = vt.R();
                    scalar sumDiag = tr(R);
                    scalar sumMagDiag = mag(R.xx())+mag(R.yy())+mag(R.zz());

                    if (mag(sumMagDiag-3) < tol && mag(sumDiag+1) < tol)
                    {
                        antiCyclic = true;
                    }
                }

                if (antiCyclic)
                {
                    // 180 degree rotational. Reset transformation.
                    permutation[matchTransI] = 0;
                }
                else
                {
                    FatalErrorIn
                    (
                        "Foam::label "
                        "Foam::globalIndexAndTransform::addToTransformIndex\n"
                        "(\n"
                            "const label,\n"
                            "const label,\n"
                            "const bool\n"
                        ") const\n"
                    )   << "More than one patch accessing the same transform "
                        << "but not of the same sign." << endl
                        << "patch:" << mesh_.boundaryMesh()[patchI].name()
                        << " transform:" << matchTransI << " sign:" << sign
                        << "  current transforms:" << permutation
                        << exit(FatalError);
                }
            }
            else
            {
                permutation[matchTransI] = 0;
            }
        }
        else
        {
            permutation[matchTransI] = sign;
        }


        // Re-encode permutation
        // ~~~~~~~~~~~~~~~~~~~~~

        return encodeTransformIndex(permutation);
    }
    else
    {
        return transformIndex;
    }
}


Foam::label Foam::globalIndexAndTransform::minimumTransformIndex
(
    const label transformIndex0,
    const label transformIndex1
) const
{
    if (transformIndex0 == transformIndex1)
    {
        return transformIndex0;
    }


    // Count number of transforms
    FixedList<label, 3> permutation0 = decodeTransformIndex(transformIndex0);
    label n0 = 0;
    forAll(permutation0, i)
    {
        if (permutation0[i] != 0)
        {
            n0++;
        }
    }

    FixedList<label, 3> permutation1 = decodeTransformIndex(transformIndex1);
    label n1 = 0;
    forAll(permutation1, i)
    {
        if (permutation1[i] != 0)
        {
            n1++;
        }
    }

    if (n0 <= n1)
    {
        return transformIndex0;
    }
    else
    {
        return transformIndex1;
    }
}


Foam::label Foam::globalIndexAndTransform::subtractTransformIndex
(
    const label transformIndex0,
    const label transformIndex1
) const
{
    FixedList<label, 3> permutation0 = decodeTransformIndex(transformIndex0);
    FixedList<label, 3> permutation1 = decodeTransformIndex(transformIndex1);

    forAll(permutation0, i)
    {
        permutation0[i] -= permutation1[i];
    }

    return encodeTransformIndex(permutation0);
}


Foam::labelPair Foam::globalIndexAndTransform::encode
(
    const label index,
    const label transformIndex
)
{
    return encode(Pstream::myProcNo(), index, transformIndex);
}


Foam::labelPair Foam::globalIndexAndTransform::encode
(
    const label procI,
    const label index,
    const label transformIndex
)
{
    if (transformIndex < 0 || transformIndex >= base_)
    {
        FatalErrorIn
        (
            "Foam::labelPair Foam::globalIndexAndTransform::encode"
            "("
                "const label procI, "
                "const label index, "
                "const label transformIndex"
            ")"
        )
            << "TransformIndex " << transformIndex
            << " is outside allowed range of 0 to "
            << base_ - 1
            << abort(FatalError);
    }

    if (procI > labelMax/base_)
    {
        FatalErrorIn
        (
            "Foam::labelPair Foam::globalIndexAndTransform::encode"
            "("
                "const label procI, "
                "const label index, "
                "const label transformIndex"
            ")"
        )
            << "Overflow : encoding processor " << procI << " in base " << base_
            << " exceeds capability of label (" << labelMax
            << "). Please recompile with larger datatype for label."
            << exit(FatalError);
    }

    return labelPair
    (
        index,
        transformIndex + procI*base_
    );
}


Foam::label Foam::globalIndexAndTransform::index
(
    const labelPair& globalIAndTransform
)
{
    return globalIAndTransform.first();
}


Foam::label Foam::globalIndexAndTransform::processor
(
    const labelPair& globalIAndTransform
)
{
    return globalIAndTransform.second()/base_;
}


Foam::label Foam::globalIndexAndTransform::transformIndex
(
    const labelPair& globalIAndTransform
)
{
    return globalIAndTransform.second() % base_;
}


Foam::label Foam::globalIndexAndTransform::nIndependentTransforms() const
{
    return transforms_.size();
}


const Foam::List<Foam::vectorTensorTransform>&
Foam::globalIndexAndTransform::transforms() const
{
    return transforms_;
}


const Foam::List<Foam::vectorTensorTransform>&
Foam::globalIndexAndTransform::transformPermutations() const
{
    return transformPermutations_;
}


Foam::label Foam::globalIndexAndTransform::nullTransformIndex() const
{
    return nullTransformIndex_;
}


const Foam::List<Foam::Pair<Foam::label> >&
Foam::globalIndexAndTransform::patchTransformSign() const
{
    return patchTransformSign_;
}


const Foam::vectorTensorTransform& Foam::globalIndexAndTransform::transform
(
    label transformIndex
) const
{
    return transformPermutations_[transformIndex];
}


Foam::labelList Foam::globalIndexAndTransform::transformIndicesForPatches
(
    const labelHashSet& patchIs
) const
{
    List<label> permutation(transforms_.size(), 0);

    labelList selectedTransformIs(0);

    if (patchIs.empty() || transforms_.empty())
    {
        return selectedTransformIs;
    }

    forAllConstIter(labelHashSet, patchIs, iter)
    {
        label patchI = iter.key();

        const Pair<label>& transSign = patchTransformSign_[patchI];

        label matchTransI = transSign.first();

        if (matchTransI > -1)
        {
            label sign = transSign.second();

            // If this transform been found already by a patch?
            if (permutation[matchTransI] != 0)
            {
                // If so, if they have opposite signs, then this is
                // considered an error.  They are allowed to be the
                // same sign, but this only results in a single
                // transform.
                if (permutation[matchTransI] != sign)
                {
                    FatalErrorIn
                    (
                        "const Foam::List<Foam::vectorTensorTransform>& "
                        "Foam::globalIndexAndTransform::transformsForPatches"
                        "("
                            "const labelList& patchIs"
                        ") const"
                    )
                        << "More than one patch accessing the same transform "
                        << "but not of the same sign."
                        << exit(FatalError);
                }
            }
            else
            {
                permutation[matchTransI] = sign;
            }
        }
    }

    label nUsedTrans = round(sum(mag(permutation)));

    if (nUsedTrans == 0)
    {
        return selectedTransformIs;
    }

    // Number of selected transformations
    label nSelTrans = pow(label(2), nUsedTrans) - 1;

    // Pout<< nl << permutation << nl << endl;

    selectedTransformIs.setSize(nSelTrans);

    switch (nUsedTrans)
    {
        case 1:
        {
            selectedTransformIs[0] = encodeTransformIndex(permutation);

            break;
        }
        case 2:
        {
            List<label> tempPermutation = permutation;

            label a = 0;
            label b = 1;

            // When there are two selected transforms out of three, we
            // need to choose which of them are being permuted
            if (transforms_.size() > nUsedTrans)
            {
                if (permutation[0] == 0)
                {
                    a = 1;
                    b = 2;
                }
                else if (permutation[1] == 0)
                {
                    a = 0;
                    b = 2;
                }
                else if (permutation[2] == 0)
                {
                    a = 0;
                    b = 1;
                }
            }

            tempPermutation[a] = a;
            tempPermutation[b] = permutation[b];

            selectedTransformIs[0] = encodeTransformIndex(tempPermutation);

            tempPermutation[a] = permutation[a];
            tempPermutation[b] = a;

            selectedTransformIs[1] = encodeTransformIndex(tempPermutation);

            tempPermutation[a] = permutation[a];
            tempPermutation[b] = permutation[b];

            selectedTransformIs[2] = encodeTransformIndex(tempPermutation);

            break;
        }
        case 3:
        {
            List<label> tempPermutation = permutation;

            tempPermutation[0] = 0;
            tempPermutation[1] = 0;
            tempPermutation[2] = permutation[2];

            selectedTransformIs[0] = encodeTransformIndex(tempPermutation);

            tempPermutation[0] = 0;
            tempPermutation[1] = permutation[1];
            tempPermutation[2] = 0;

            selectedTransformIs[1] = encodeTransformIndex(tempPermutation);

            tempPermutation[0] = 0;
            tempPermutation[1] = permutation[1];
            tempPermutation[2] = permutation[2];

            selectedTransformIs[2] = encodeTransformIndex(tempPermutation);

            tempPermutation[0] = permutation[0];
            tempPermutation[1] = 0;
            tempPermutation[2] = 0;

            selectedTransformIs[3] = encodeTransformIndex(tempPermutation);

            tempPermutation[0] = permutation[0];
            tempPermutation[1] = 0;
            tempPermutation[2] = permutation[2];

            selectedTransformIs[4] = encodeTransformIndex(tempPermutation);

            tempPermutation[0] = permutation[0];
            tempPermutation[1] = permutation[1];
            tempPermutation[2] = 0;

            selectedTransformIs[5] = encodeTransformIndex(tempPermutation);

            tempPermutation[0] = permutation[0];
            tempPermutation[1] = permutation[1];
            tempPermutation[2] = permutation[2];

            selectedTransformIs[6] = encodeTransformIndex(tempPermutation);

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "const Foam::List<Foam::vectorTensorTransform>& "
                "Foam::globalIndexAndTransform::transformsForPatches"
                "("
                    "const labelList& patchIs"
                ") const"
            )
                << "Only 1-3 transforms are possible."
                << exit(FatalError);
        }
    }

    return selectedTransformIs;
}


Foam::pointField Foam::globalIndexAndTransform::transformPatches
(
    const labelHashSet& patchIs,
    const point& pt
) const
{
    labelList transIs = transformIndicesForPatches(patchIs);

    // Pout<< patchIs << nl << transIs << endl;

    pointField transPts(transIs.size());

    forAll(transIs, tII)
    {
        transPts[tII] = transformPermutations_[transIs[tII]].transformPosition
        (
            pt
        );
    }

    return transPts;
}


// ************************************************************************* //
