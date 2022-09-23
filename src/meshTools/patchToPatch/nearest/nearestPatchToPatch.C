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
    defineTypeNameAndDebug(nearest, 0);
    addToRunTimeSelectionTable(patchToPatch, nearest, bool);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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
    if
    (
        nearby::intersectFaces
        (
            srcPatch,
            srcPointNormals,
            srcPointNormals0,
            tgtPatch,
            srcFacei,
            tgtFacei
        )
    )
    {
        const scalar dSqr =
            magSqr
            (
                srcPatch.faceCentres()[srcFacei]
              - tgtPatch.faceCentres()[tgtFacei]
            );

        if (dSqr < srcDistances_[srcFacei])
        {
            srcDistances_[srcFacei] = dSqr;
            Swap
            (
                srcLocalTgtFaces_[srcFacei].first(),
                srcLocalTgtFaces_[srcFacei].last()
            );
        }

        if (dSqr < tgtDistances_[tgtFacei])
        {
            tgtDistances_[tgtFacei] = dSqr;
            Swap
            (
                tgtLocalSrcFaces_[tgtFacei].first(),
                tgtLocalSrcFaces_[tgtFacei].last()
            );
        }

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::patchToPatches::nearest::initialise
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

    srcDistances_.resize(srcPatch.size());
    srcDistances_ = vGreat;

    tgtDistances_.resize(tgtPatch.size());
    tgtDistances_ = vGreat;
}


Foam::labelList Foam::patchToPatches::nearest::finaliseLocal
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch
)
{
    const labelList newToOldLocalTgtFace =
        nearby::finaliseLocal
        (
            srcPatch,
            srcPointNormals,
            srcPointNormals0,
            tgtPatch
        );

    tgtDistances_ = List<scalar>(tgtDistances_, newToOldLocalTgtFace);

    return newToOldLocalTgtFace;
}


void Foam::patchToPatches::nearest::rDistributeTgt
(
    const primitiveOldTimePatch& tgtPatch
)
{
    // Keep only the closest opposing face
    forAll(srcLocalTgtFaces_, srcFacei)
    {
        srcLocalTgtFaces_[srcFacei].resize
        (
            min(srcLocalTgtFaces_[srcFacei].size(), 1)
        );
    }
    forAll(tgtLocalSrcFaces_, tgtFacei)
    {
        tgtLocalSrcFaces_[tgtFacei].resize
        (
            min(tgtLocalSrcFaces_[tgtFacei].size(), 1)
        );
    }

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
    nearby::rDistributeTgt(tgtPatch);

    // Reverse distribute the distances
    patchToPatchTools::rDistributeListList
    (
        tgtPatch.size(),
        tgtMapPtr_(),
        tgtDistances
    );

    // If there is more than one address, remove all but the closest
    tgtDistances_.resize(tgtLocalSrcFaces_.size());
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


Foam::label Foam::patchToPatches::nearest::finalise
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch,
    const transformer& tgtToSrc
)
{
    // Keep only the closest opposing face
    forAll(srcLocalTgtFaces_, srcFacei)
    {
        srcLocalTgtFaces_[srcFacei].resize
        (
            min(srcLocalTgtFaces_[srcFacei].size(), 1)
        );
    }
    forAll(tgtLocalSrcFaces_, tgtFacei)
    {
        tgtLocalSrcFaces_[tgtFacei].resize
        (
            min(tgtLocalSrcFaces_[tgtFacei].size(), 1)
        );
    }

    const label nCouples =
        nearby::finalise
        (
            srcPatch,
            srcPointNormals,
            srcPointNormals0,
            tgtPatch,
            tgtToSrc
        );

    if (debug)
    {
        auto countNonEmpty = [](const List<DynamicList<label>>& ll)
        {
            label result = 0;

            forAll(ll, i)
            {
                result += !ll[i].empty();
            }

            return returnReduce(result, sumOp<label>());
        };

        Info<< indent
            << "Coupled " << countNonEmpty(srcLocalTgtFaces_)
            << "/" << returnReduce(srcLocalTgtFaces_.size(), sumOp<label>())
            << " source faces and " << countNonEmpty(tgtLocalSrcFaces_)
            << "/" << returnReduce(tgtLocalSrcFaces_.size(), sumOp<label>())
            << " target faces" << endl;
    }

    if (debug && !Pstream::parRun())
    {
        auto writeConnections = []
        (
            const primitivePatch& patch,
            const primitivePatch& otherPatch,
            const bool isSrc,
            const List<DynamicList<label>>& faceLocalOtherFaces
        )
        {
            const word name =
                typeName + "_" + (isSrc ? "src" : "tgt") + "Connections";

            OBJstream obj(name + ".obj");

            forAll(faceLocalOtherFaces, facei)
            {
                const point& p = patch.faceCentres()[facei];
                forAll(faceLocalOtherFaces[facei], i)
                {
                    const label otherFacei = faceLocalOtherFaces[facei][i];
                    const point& q = otherPatch.faceCentres()[otherFacei];
                    obj.write(linePointRef(p, q));
                }
            }
        };

        writeConnections(srcPatch, tgtPatch, true, srcLocalTgtFaces_);
        writeConnections(tgtPatch, srcPatch, false, tgtLocalSrcFaces_);

        auto writeNotConnected = []
        (
            const primitivePatch& patch,
            const List<DynamicList<label>>& faceLocalOtherFaces,
            const bool isSrc
        )
        {
            DynamicList<label> unconnected;
            forAll(faceLocalOtherFaces, facei)
            {
                if (faceLocalOtherFaces[facei].empty())
                {
                    unconnected.append(facei);
                }
            }

            const word name =
                typeName + "_" + (isSrc ? "src" : "tgt") + "NotConnected";

            vtkWritePolyData::write
            (
                name + ".vtk",
                name,
                false,
                patch.localPoints(),
                labelList(),
                labelListList(),
                UIndirectList<face>(patch.localFaces(), unconnected)
            );
        };

        writeNotConnected(srcPatch, srcLocalTgtFaces_, true);
        writeNotConnected(tgtPatch, tgtLocalSrcFaces_, false);
    }

    return nCouples;
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
    nearby(reverse),
    srcDistances_(),
    tgtDistances_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToPatches::nearest::~nearest()
{}


// ************************************************************************* //
