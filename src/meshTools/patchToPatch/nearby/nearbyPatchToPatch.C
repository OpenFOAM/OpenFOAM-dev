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

#include "nearbyPatchToPatch.H"
#include "boundSphere.H"
#include "OFstream.H"
#include "OBJstream.H"
#include "vtkWritePolyData.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatches
{
    defineTypeNameAndDebug(nearby, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::treeBoundBox Foam::patchToPatches::nearby::srcBox
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


bool Foam::patchToPatches::nearby::intersectFaces
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch,
    const label srcFacei,
    const label tgtFacei
)
{
    const point& srcC = srcSpheres_[srcFacei].first();
    const scalar srcR = srcSpheres_[srcFacei].second();
    const point& tgtC = tgtSpheres_[tgtFacei].first();
    const scalar tgtR = tgtSpheres_[tgtFacei].second();

    if (magSqr(srcC - tgtC) < sqr(srcR + tgtR))
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


void Foam::patchToPatches::nearby::initialise
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

    srcSpheres_.resize(srcPatch.size());
    forAll(srcPatch, srcFacei)
    {
        srcSpheres_[srcFacei] =
            boundSphere
            (
                UIndirectList<point>(srcPatch.points(), srcPatch[srcFacei])
            );
    }

    tgtSpheres_.resize(tgtPatch.size());
    forAll(tgtPatch, tgtFacei)
    {
        tgtSpheres_[tgtFacei] =
            boundSphere
            (
                UIndirectList<point>(tgtPatch.points(), tgtPatch[tgtFacei])
            );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatches::nearby::nearby(const bool reverse)
:
    patchToPatch(reverse),
    srcSpheres_(),
    tgtSpheres_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToPatches::nearby::~nearby()
{}


// ************************************************************************* //
