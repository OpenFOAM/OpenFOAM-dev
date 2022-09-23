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


// * * * * * * * * * * * Private Static Member Functions * * * * * * * * * * //

bool Foam::patchToPatches::inverseDistance::rayHitsFace
(
    const point& p,
    const vector& r,
    const face& f,
    const pointField& ps
)
{
    using namespace constant::mathematical;

    const tensor T = tensor::I - sqr(r);

    scalar angle = 0;

    forAll(f, i)
    {
        const vector& a = T & (ps[f[i]] - p);
        const vector& b = T & (ps[f[f.fcIndex(i)]] - p);

        const scalar meanMagSqrAB = (magSqr(a) + magSqr(b))/2;
        const scalar geometricMeanMagSqrAB = sqrt(magSqr(a)*magSqr(b));

        // This indicates that we have hit a point to within round off error
        if (geometricMeanMagSqrAB < small*meanMagSqrAB) return true;

        angle -=
            sign(r & (a ^ b))
           *acos(min(max(-1, (a & b)/geometricMeanMagSqrAB), +1));
    }

    return pi < angle && angle < 3*pi;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchToPatches::inverseDistance::initialise
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch
)
{
    nearby::initialise
    (
        srcPatch,
        srcPointNormals,
        srcPointNormals0,
        tgtPatch
    );

    srcWeights_.resize(srcPatch.size());
    forAll(srcWeights_, i)
    {
        srcWeights_[i].clear();
    }

    tgtWeights_.resize(tgtPatch.size());
    forAll(tgtWeights_, i)
    {
        tgtWeights_[i].clear();
    }
}


void Foam::patchToPatches::inverseDistance::generateWeights
(
    const primitiveOldTimePatch& srcPatch,
    const primitiveOldTimePatch& tgtPatch
)
{
    auto generate = []
    (
        const primitiveOldTimePatch& patch,
        const primitiveOldTimePatch& otherPatch,
        const bool reverse,
        List<DynamicList<label>>& otherFaces,
        List<DynamicList<scalar>>& weights
    )
    {
        forAll(otherFaces, facei)
        {
            if (otherFaces[facei].empty()) continue;

            label otherFacei = -1;

            // Find the other face that "contains" this face's centre
            forAll(otherFaces[facei], i)
            {
                if
                (
                    rayHitsFace
                    (
                        patch.faceCentres()[facei],
                        (reverse ? -1 : +1)*patch.faceNormals()[facei],
                        otherPatch[otherFaces[facei][i]],
                        otherPatch.points()
                    )
                )
                {
                    otherFacei = otherFaces[facei][i];
                    break;
                }
            }

            const point& c = patch.faceCentres()[facei];

            // If the above failed, find the closest
            if (otherFacei == -1)
            {
                scalar minDistSqr = vGreat;

                forAll(otherFaces[facei], i)
                {
                    const point& otherC =
                        otherPatch.faceCentres()[otherFaces[facei][i]];
                    const scalar distSqr = magSqr(c - otherC);
                    if (distSqr < minDistSqr)
                    {
                        minDistSqr = distSqr;
                        otherFacei = otherFaces[facei][i];
                    }
                }
            }

            // Remove all faces
            otherFaces[facei].clear();

            // Add the found face and all its neighbours
            const point& otherC = otherPatch.faceCentres()[otherFacei];
            otherFaces[facei].append(otherFacei);
            weights[facei].append(1/(mag(c - otherC) + rootVSmall));

            forAll(otherPatch.faceFaces()[otherFacei], i)
            {
                const label otherFacej = otherPatch.faceFaces()[otherFacei][i];

                const point& otherC = otherPatch.faceCentres()[otherFacej];
                otherFaces[facei].append(otherFacej);
                weights[facei].append(1/(mag(c - otherC) + rootVSmall));
            }
        }
    };

    generate(srcPatch, tgtPatch, reverse_, srcLocalTgtFaces_, srcWeights_);
    generate(tgtPatch, srcPatch, reverse_, tgtLocalSrcFaces_, tgtWeights_);
}


Foam::labelList Foam::patchToPatches::inverseDistance::finaliseLocal
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch
)
{
    // Transfer weight addressing into actual addressing
    generateWeights(srcPatch, tgtPatch);

    const labelList newToOldLocalTgtFace =
        nearby::finaliseLocal
        (
            srcPatch,
            srcPointNormals,
            srcPointNormals0,
            tgtPatch
        );

    tgtWeights_ = List<DynamicList<scalar>>(tgtWeights_, newToOldLocalTgtFace);

    return newToOldLocalTgtFace;
}


void Foam::patchToPatches::inverseDistance::rDistributeTgt
(
    const primitiveOldTimePatch& tgtPatch
)
{
    // Let the base class reverse distribute the addressing
    nearby::rDistributeTgt(tgtPatch);

    // Reverse distribute the weights
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

    // Transfer weight addressing into actual addressing (if not done in the
    // finaliseLocal method above)
    if (isSingleProcess())
    {
        generateWeights(srcPatch, tgtPatch);
    }

    // Normalise the weights
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

        Info<< indent
            << "Number of source faces by number of target connections = "
            << histogram(srcLocalTgtFaces_) << nl
            << indent
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
    nearby(reverse),
    srcWeights_(),
    tgtWeights_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToPatches::inverseDistance::~inverseDistance()
{}


// ************************************************************************* //
