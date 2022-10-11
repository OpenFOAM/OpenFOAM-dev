/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "inverseDistancePatchToPatch.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatches
{
    defineTypeNameAndDebug(inverseDistance, 0);
    addToRunTimeSelectionTable(patchToPatch, inverseDistance, bool);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::treeBoundBox Foam::patchToPatches::inverseDistance::srcBox
(
    const face& srcFace,
    const pointField& srcPoints,
    const vectorField& srcPointNormals
) const
{
    const treeBoundBox bb(srcPoints, srcFace);

    const point c = bb.midpoint();
    const scalar l = bb.maxDim();

    return treeBoundBox(c - l*vector::one, c + l*vector::one);
}


bool Foam::patchToPatches::inverseDistance::inside
(
    const face& f,
    const pointField& ps,
    const point& p,
    const vector& r
) const
{
    using namespace constant::mathematical;

    const tensor T = tensor::I - sqr(r);

    scalar angle = 0;

    forAll(f, i)
    {
        const vector& a = T & (ps[f[i]] - p);
        const vector& b = T & (ps[f[f.fcIndex(i)]] - p);
        const scalar magAB = sqrt(magSqr(a)*magSqr(b));
        angle +=
            (reverse() ? +1 : -1)
           *sign(r & (a ^ b))
           *acos(min(max(-1, (a & b)/magAB), +1));
    }

    return pi < angle && angle < 3*pi;
}


bool Foam::patchToPatches::inverseDistance::intersectFaces
(
    const primitivePatch& patch,
    const primitivePatch& otherPatch,
    const label facei,
    const label otherFacei,
    DynamicList<label>& faceOtherFaces,
    DynamicList<scalar>& faceWeights
) const
{
    const face& f = otherPatch[otherFacei];
    const pointField& ps = otherPatch.points();

    const point& p = patch.faceCentres()[facei];
    const vector& r = patch.faceNormals()[facei];

    bool intersectsOther = inside(f, ps, p, r);

    if (!intersectsOther)
    {
        forAll(otherPatch.faceFaces()[otherFacei], otherFaceFacei)
        {
            const label otherFacej =
                otherPatch.faceFaces()[otherFacei][otherFaceFacei];

            const face& g = otherPatch[otherFacej];

            if (inside(g, ps, p, r))
            {
                intersectsOther = true;
                break;
            }
        }
    }

    if (intersectsOther)
    {
        faceOtherFaces.append(otherFacei);
        faceWeights.append
        (
            1/max(mag(p - otherPatch.faceCentres()[otherFacei]), vSmall)
        );
    }

    return intersectsOther;
}


bool Foam::patchToPatches::inverseDistance::intersectFaces
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch,
    const label srcFacei,
    const label tgtFacei
)
{
    const bool srcCouples =
        intersectFaces
        (
            srcPatch,
            tgtPatch,
            srcFacei,
            tgtFacei,
            srcLocalTgtFaces_[srcFacei],
            srcWeights_[srcFacei]
        );

    const bool tgtCouples =
        intersectFaces
        (
            tgtPatch,
            srcPatch,
            tgtFacei,
            srcFacei,
            tgtLocalSrcFaces_[tgtFacei],
            tgtWeights_[tgtFacei]
        );

    return srcCouples || tgtCouples;
}


void Foam::patchToPatches::inverseDistance::initialise
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch
)
{
    patchToPatch::initialise
    (
        srcPatch,
        srcPointNormals,
        srcPointNormals0,
        tgtPatch
    );

    srcWeights_.resize(srcPatch.size());
    tgtWeights_.resize(tgtPatch.size());
}


Foam::labelList Foam::patchToPatches::inverseDistance::trimLocalTgt
(
    const primitiveOldTimePatch& localTgtPatch
)
{
    const labelList newToOldLocalTgtFace =
        patchToPatch::trimLocalTgt(localTgtPatch);

    tgtWeights_ = List<DynamicList<scalar>>(tgtWeights_, newToOldLocalTgtFace);

    return newToOldLocalTgtFace;
}


void Foam::patchToPatches::inverseDistance::rDistributeTgt
(
    const primitiveOldTimePatch& tgtPatch
)
{
    patchToPatch::rDistributeTgt(tgtPatch);

    patchToPatchTools::rDistributeListList
    (
        tgtPatch.size(),
        tgtMapPtr_(),
        tgtWeights_
    );
}


Foam::label Foam::patchToPatches::inverseDistance::finalise
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch,
    const transformer& tgtToSrc
)
{
    const label nCouples =
        patchToPatch::finalise
        (
            srcPatch,
            srcPointNormals,
            srcPointNormals0,
            tgtPatch,
            tgtToSrc
        );

    forAll(srcWeights_, srcFacei)
    {
        const scalar w = sum(srcWeights_[srcFacei]);

        forAll(srcWeights_[srcFacei], i)
        {
            srcWeights_[srcFacei][i] /= max(w, vSmall);
        }
    }

    forAll(tgtWeights_, tgtFacei)
    {
        const scalar w = sum(tgtWeights_[tgtFacei]);

        forAll(tgtWeights_[tgtFacei], i)
        {
            tgtWeights_[tgtFacei][i] /= max(w, vSmall);
        }
    }

    if (debug)
    {
        auto histogram = [](const List<DynamicList<label>>& ll)
        {
            labelList result;
            forAll(ll, i)
            {
                result.setSize(max(result.size(), ll[i].size() + 1), 0);
                result[ll[i].size()] ++;
            }

            result.resize(returnReduce(result.size(), maxOp<label>()), 0);

            Pstream::listCombineGather(result, plusEqOp<label>());
            Pstream::listCombineScatter(result);

            return result;
        };

        Info<< "Number of source faces by number of target connections = "
            << histogram(srcLocalTgtFaces_) << nl
            << "Number of target faces by number of source connections = "
            << histogram(tgtLocalSrcFaces_) << endl;
    }

    return nCouples;
}


Foam::tmpNrc<Foam::List<Foam::DynamicList<Foam::scalar>>>
Foam::patchToPatches::inverseDistance::srcWeights() const
{
    return srcWeights_;
}


Foam::tmpNrc<Foam::List<Foam::DynamicList<Foam::scalar>>>
Foam::patchToPatches::inverseDistance::tgtWeights() const
{
    return tgtWeights_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatches::inverseDistance::inverseDistance(const bool reverse)
:
    patchToPatch(reverse),
    srcWeights_(),
    tgtWeights_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToPatches::inverseDistance::~inverseDistance()
{}


// ************************************************************************* //
