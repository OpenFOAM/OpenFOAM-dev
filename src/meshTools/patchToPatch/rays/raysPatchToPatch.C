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

#include "raysPatchToPatch.H"
#include "intersectionPatchToPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatches
{
    defineTypeNameAndDebug(rays, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::treeBoundBox Foam::patchToPatches::rays::srcBox
(
    const face& srcFace,
    const pointField& srcPoints,
    const vectorField& srcPointNormals
) const
{
    return intersection::srcBoxStatic(srcFace, srcPoints, srcPointNormals);
}


bool Foam::patchToPatches::rays::intersectFaces
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch,
    const label srcFacei,
    const label tgtFacei
)
{
    const treeBoundBox srcFaceBox =
        patchToPatch::srcBox
        (
            srcPatch,
            srcPointNormals,
            srcPointNormals0,
            srcFacei
        );
    const treeBoundBox tgtFaceBox =
        patchToPatch::tgtBox(tgtPatch, tgtFacei);

    if (srcFaceBox.overlaps(tgtFaceBox))
    {
        srcLocalTgtFaces_[srcFacei].append(tgtFacei);
        tgtLocalSrcFaces_[tgtFacei].append(srcFacei);

        return true;
    }
    else
    {
        return false;
    }
}


Foam::labelList Foam::patchToPatches::rays::finaliseLocal
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch
)
{
    const labelList newToOldLocalTgtFace =
        patchToPatch::finaliseLocal
        (
            srcPatch,
            srcPointNormals,
            srcPointNormals0,
            tgtPatch
        );

    localTgtPatchPtr_.reset
    (
        new PrimitiveOldTimePatch<faceList, pointField>
        (
            faceList(tgtPatch, newToOldLocalTgtFace),
            tgtPatch.points(),
            tgtPatch.points0()
        )
    );

    return newToOldLocalTgtFace;
}


void Foam::patchToPatches::rays::distributeSrc
(
    const primitiveOldTimePatch& srcPatch
)
{
    localSrcProcFacesPtr_.reset
    (
        new List<remote>
        (
            distributePatch(srcMapPtr_(), srcPatch, localSrcPatchPtr_)
        )
    );
}


Foam::label Foam::patchToPatches::rays::finalise
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

    // Transform the source-local target patch back to the target
    if (!isSingleProcess() && !isNull(tgtToSrc))
    {
        autoPtr<PrimitiveOldTimePatch<faceList, pointField>>
            localTgtPatchPtr(localTgtPatchPtr_.ptr());

        localTgtPatchPtr_.set
        (
            new PrimitiveOldTimePatch<faceList, pointField>
            (
                localTgtPatchPtr(),
                tgtToSrc.invTransformPosition(localTgtPatchPtr().points()),
                isNull(localTgtPatchPtr().points0())
              ? NullObjectRef<pointField>()
              : tgtToSrc.invTransformPosition(localTgtPatchPtr().points0())()
            )
        );
    }

    return nCouples;
}


Foam::remote Foam::patchToPatches::rays::ray
(
    const primitiveOldTimePatch& outPatch,
    const autoPtr<PrimitiveOldTimePatch<faceList, pointField>>&
        localOutPatchPtr,
    const autoPtr<List<remote>>& localOutProcFacesPtr,
    const List<DynamicList<label>>& inLocalOutFaces,
    const scalar fraction,
    const label inFacei,
    const point& inP,
    const vector& inN,
    point& outP
) const
{
    forAll(inLocalOutFaces[inFacei], i)
    {
        const label localOutFacei = inLocalOutFaces[inFacei][i];

        const face& outF =
            isSingleProcess()
          ? outPatch[localOutFacei]
          : localOutPatchPtr()[localOutFacei];

        const pointField& outPoints =
            isSingleProcess()
          ? outPatch.points()
          : localOutPatchPtr().points();
        const pointField& outPoints0 =
            isSingleProcess()
          ? outPatch.points0()
          : localOutPatchPtr().points0();

        const pointField outPoly
        (
            (1 - fraction)*outF.points(outPoints0)
          + fraction*outF.points(outPoints)
        );

        const pointHit ray =
            face(identity(outPoly.size()))
           .ray(inP, inN, outPoly, Foam::intersection::algorithm::visible);

        if (ray.hit())
        {
            outP = ray.rawPoint();

            return
                isSingleProcess()
              ? remote({Pstream::myProcNo(), localOutFacei})
              : localOutProcFacesPtr()[localOutFacei];
        }
    }

    return remote({-1, -1});
}


Foam::tmpNrc<Foam::List<Foam::DynamicList<Foam::scalar>>>
Foam::patchToPatches::rays::srcWeights() const
{
    NotImplemented;
    return tmpNrc<List<DynamicList<scalar>>>(nullptr);
}


Foam::tmpNrc<Foam::List<Foam::DynamicList<Foam::scalar>>>
Foam::patchToPatches::rays::tgtWeights() const
{
    NotImplemented;
    return tmpNrc<List<DynamicList<scalar>>>(nullptr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatches::rays::rays
(
    const bool reverse
)
:
    patchToPatch(reverse)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToPatches::rays::~rays()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::remote Foam::patchToPatches::rays::srcToTgtRay
(
    const primitiveOldTimePatch& tgtPatch,
    const scalar fraction,
    const label srcFacei,
    const vector& srcP,
    const vector& srcN,
    point& tgtP
) const
{
    return
        ray
        (
            tgtPatch,
            localTgtPatchPtr_,
            localTgtProcFacesPtr_,
            srcLocalTgtFaces_,
            fraction,
            srcFacei,
            srcP,
            srcN,
            tgtP
        );
}


Foam::remote Foam::patchToPatches::rays::tgtToSrcRay
(
    const primitiveOldTimePatch& srcPatch,
    const scalar fraction,
    const label tgtFacei,
    const vector& tgtP,
    const vector& tgtN,
    point& srcP
) const
{
    return
        ray
        (
            srcPatch,
            localSrcPatchPtr_,
            localSrcProcFacesPtr_,
            tgtLocalSrcFaces_,
            fraction,
            tgtFacei,
            tgtP,
            tgtN,
            srcP
        );
}


// ************************************************************************* //
