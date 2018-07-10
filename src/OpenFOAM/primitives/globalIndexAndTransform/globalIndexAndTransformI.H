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

#include "polyMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::globalIndexAndTransform::less::operator()
(
    const labelPair& a,
    const labelPair& b
) const
{
    label procA = gi_.processor(a);
    label procB = gi_.processor(b);

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
        label indexA = gi_.index(a);
        label indexB = gi_.index(b);

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
            label transformA = gi_.transformIndex(a);
            label transformB = gi_.transformIndex(b);

            return transformA < transformB;
        }
    }
}


Foam::label Foam::globalIndexAndTransform::encodeTransformIndex
(
    const labelList& permutationIndices
) const
{
    if (permutationIndices.size() != transforms_.size())
    {
        FatalErrorInFunction
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
            FatalErrorInFunction
                << "permutationIndices " << permutationIndices
                << "are illegal, they must all be only -1, 0 or +1"
                << abort(FatalError);
        }

        transformIndex += (permutationIndices[b] + 1)*w;

        w *= 3;
    }

    return transformIndex;
}


Foam::labelList Foam::globalIndexAndTransform::decodeTransformIndex
(
    const label transformIndex
) const
{
    labelList permutation(transforms_.size(), 0);

    label t = transformIndex;
    forAll(permutation, i)
    {
        permutation[i] = (t%3)-1;
        t /= 3;
    }

    return permutation;
}


Foam::label Foam::globalIndexAndTransform::addToTransformIndex
(
    const label transformIndex,
    const label patchi,
    const bool isSendingSide,
    const scalar tol
) const
{
    const labelPair& transSign = patchTransformSign_[patchi];

    label matchTransI = transSign.first();

    if (matchTransI >= transforms_.size())
    {
        FatalErrorInFunction
            << "patch:" << mesh_.boundaryMesh()[patchi].name()
            << " transform:" << matchTransI
            << " out of possible transforms:" << transforms_
            << exit(FatalError);
        return labelMin;
    }
    else if (matchTransI == -1)
    {
        // No additional transformation for this patch
        return transformIndex;
    }
    else
    {
        // Decode current set of transforms
        labelList permutation(decodeTransformIndex(transformIndex));


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
                FatalErrorInFunction
                    << "patch:" << mesh_.boundaryMesh()[patchi].name()
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
                if (mag(vt.t()) < small && vt.hasR())
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
                    FatalErrorInFunction
                        << "More than one patch accessing the same transform "
                        << "but not of the same sign." << endl
                        << "patch:" << mesh_.boundaryMesh()[patchi].name()
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
    labelList permutation0(decodeTransformIndex(transformIndex0));
    label n0 = 0;
    forAll(permutation0, i)
    {
        if (permutation0[i] != 0)
        {
            n0++;
        }
    }

    labelList permutation1(decodeTransformIndex(transformIndex1));
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
    labelList permutation0(decodeTransformIndex(transformIndex0));
    labelList permutation1(decodeTransformIndex(transformIndex1));

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
) const
{
    return encode(Pstream::myProcNo(), index, transformIndex);
}


Foam::labelPair Foam::globalIndexAndTransform::encode
(
    const label proci,
    const label index,
    const label transformIndex
) const
{
    if (transformIndex < 0 || transformIndex >= transformPermutations_.size())
    {
        FatalErrorInFunction
            << "TransformIndex " << transformIndex
            << " is outside allowed range of 0 to "
            << transformPermutations_.size() - 1
            << abort(FatalError);
    }

    if (proci > labelMax/transformPermutations_.size())
    {
        FatalErrorInFunction
            << "Overflow : encoding processor " << proci
            << " in base " << transformPermutations_.size()
            << " exceeds capability of label (" << labelMax
            << "). Please recompile with larger datatype for label."
            << exit(FatalError);
    }

    return labelPair
    (
        index,
        transformIndex + proci*transformPermutations_.size()
    );
}


Foam::label Foam::globalIndexAndTransform::index
(
    const labelPair& globalIAndTransform
) const
{
    return globalIAndTransform.first();
}


Foam::label Foam::globalIndexAndTransform::processor
(
    const labelPair& globalIAndTransform
) const
{
    return globalIAndTransform.second()/transformPermutations_.size();
}


Foam::label Foam::globalIndexAndTransform::transformIndex
(
    const labelPair& globalIAndTransform
) const
{
    return globalIAndTransform.second()%transformPermutations_.size();
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


const Foam::labelPairList&
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
    const labelHashSet& patchis
) const
{
    labelList permutation(transforms_.size(), 0);

    labelList selectedTransformIs(0);

    if (patchis.empty() || transforms_.empty())
    {
        return selectedTransformIs;
    }

    forAllConstIter(labelHashSet, patchis, iter)
    {
        label patchi = iter.key();

        const labelPair& transSign = patchTransformSign_[patchi];

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
                    FatalErrorInFunction
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
            labelList tempPermutation = permutation;

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
            labelList tempPermutation = permutation;

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
            FatalErrorInFunction
                << "Only 1-3 transforms are possible."
                << exit(FatalError);
        }
    }

    return selectedTransformIs;
}


Foam::pointField Foam::globalIndexAndTransform::transformPatches
(
    const labelHashSet& patchis,
    const point& pt
) const
{
    labelList transIs = transformIndicesForPatches(patchis);

    // Pout<< patchis << nl << transIs << endl;

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
