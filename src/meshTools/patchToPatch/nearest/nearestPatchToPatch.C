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

#include "nearestPatchToPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatches
{
    defineTypeNameAndDebug(nearest, 0);
    addToRunTimeSelectionTable(patchToPatch, nearest, bool);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::treeBoundBox Foam::patchToPatches::nearest::srcBox
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


bool Foam::patchToPatches::nearest::intersectFaces
(
    const primitivePatch& patch,
    const primitivePatch& otherPatch,
    const label facei,
    const label otherFacei,
    DynamicList<label>& faceOtherFaces,
    scalar& faceDistance
) const
{
    auto closest = [&patch,&otherPatch]
    (
        const label facei,
        const label otherFacei
    )
    {
        const point& c = patch.faceCentres()[facei];
        const point& otherC = otherPatch.faceCentres()[otherFacei];
        const scalar distSqr = magSqr(c - otherC);

        forAll(otherPatch.faceEdges()[otherFacei], otherFaceEdgei)
        {
            const label otherEdgei =
                otherPatch.faceEdges()[otherFacei][otherFaceEdgei];

            point otherNbrC;

            if (otherPatch.edgeFaces()[otherEdgei].size() == 2)
            {
                const label facej =
                    otherPatch.edgeFaces()[otherEdgei]
                    [otherPatch.edgeFaces()[otherEdgei][0] == otherFacei];

                otherNbrC = otherPatch.faceCentres()[facej];
            }
            else
            {
                const edge& e = otherPatch.edges()[otherEdgei];
                const point& p = otherPatch.localPoints()[e[0]];
                const vector dp = e.vec(otherPatch.localPoints());
                const vector n = otherPatch.faceNormals()[otherFacei] ^ dp;

                otherNbrC = p + ((tensor::I - 2*sqr(n)) & (otherC - p));
            }

            if (magSqr(c - otherNbrC) < distSqr)
            {
                return false;
            }
        }

        return true;
    };

    if (closest(facei, otherFacei))
    {
        const point& c = patch.faceCentres()[facei];
        const point& otherC = otherPatch.faceCentres()[otherFacei];
        const scalar distSqr = magSqr(c - otherC);

        if (faceOtherFaces.empty() || faceDistance > distSqr)
        {
            faceOtherFaces.clear();
            faceOtherFaces.append(otherFacei);
            faceDistance = distSqr;
        }

        return true;
    }

    const labelList& otherFaceFaces = otherPatch.faceFaces()[otherFacei];
    forAll(otherFaceFaces, otherFaceFacei)
    {
        if (closest(facei, otherFaceFaces[otherFaceFacei]))
        {
            return true;
        }
    }

    return false;
}


bool Foam::patchToPatches::nearest::intersectFaces
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
            srcDistances_[srcFacei]
        );

    const bool tgtCouples =
        intersectFaces
        (
            tgtPatch,
            srcPatch,
            tgtFacei,
            srcFacei,
            tgtLocalSrcFaces_[tgtFacei],
            tgtDistances_[tgtFacei]
        );

    return srcCouples || tgtCouples;

}


void Foam::patchToPatches::nearest::initialise
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

    srcDistances_.resize(srcPatch.size());
    srcDistances_ = vGreat;

    tgtDistances_.resize(tgtPatch.size());
    tgtDistances_ = vGreat;
}


void Foam::patchToPatches::nearest::rDistributeTgt
(
    const primitiveOldTimePatch& tgtPatch
)
{
    // Create a list-list of distances to match the addressing
    List<List<scalar>> tgtDistances(tgtLocalSrcFaces_.size());
    forAll(tgtLocalSrcFaces_, tgtFacei)
    {
        if (!tgtLocalSrcFaces_[tgtFacei].empty())
        {
            tgtDistances[tgtFacei].resize(1, tgtDistances_[tgtFacei]);
        }
    }

    // Let the base class reverse distribute the addressing
    patchToPatch::rDistributeTgt(tgtPatch);

    // Reverse distribute the distances
    rDistributeListList(tgtPatch.size(), tgtMapPtr_(), tgtDistances);

    // If there is more than one address, remove all but the closest
    forAll(tgtLocalSrcFaces_, tgtFacei)
    {
        if (tgtLocalSrcFaces_[tgtFacei].size() > 1)
        {
            const label i = findMin(tgtDistances[tgtFacei]);

            const label srcFacei = tgtLocalSrcFaces_[tgtFacei][i];

            tgtLocalSrcFaces_[tgtFacei].resize(1);
            tgtLocalSrcFaces_[tgtFacei][0] = srcFacei;

            tgtDistances_[tgtFacei] = tgtDistances[tgtFacei][i];
        }
    }
}


Foam::tmpNrc<Foam::List<Foam::DynamicList<Foam::scalar>>>
Foam::patchToPatches::nearest::srcWeights() const
{
    tmpNrc<List<DynamicList<scalar>>> tResult
    (
        new List<DynamicList<scalar>>(srcLocalTgtFaces_.size())
    );

    List<DynamicList<scalar>>& result = tResult.ref();

    forAll(srcLocalTgtFaces_, srcFacei)
    {
        if (!srcLocalTgtFaces_[srcFacei].empty())
        {
            result[srcFacei].resize(1, scalar(1));
        }
    }

    return tResult;
}


Foam::tmpNrc<Foam::List<Foam::DynamicList<Foam::scalar>>>
Foam::patchToPatches::nearest::tgtWeights() const
{
    tmpNrc<List<DynamicList<scalar>>> tResult
    (
        new List<DynamicList<scalar>>(tgtLocalSrcFaces_.size())
    );

    List<DynamicList<scalar>>& result = tResult.ref();

    forAll(tgtLocalSrcFaces_, tgtFacei)
    {
        if (!tgtLocalSrcFaces_[tgtFacei].empty())
        {
            result[tgtFacei].resize(1, scalar(1));
        }
    }

    return tResult;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatches::nearest::nearest(const bool reverse)
:
    patchToPatch(reverse)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToPatches::nearest::~nearest()
{}


// ************************************************************************* //
