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

#include "matchingPatchToPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatches
{
    defineTypeNameAndDebug(matching, 0);
    addToRunTimeSelectionTable(patchToPatch, matching, bool);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::patchToPatches::matching::finalise
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch,
    const transformer& tgtToSrc
)
{
    const label nCouples =
        nearest::finalise
        (
            srcPatch,
            srcPointNormals,
            srcPointNormals0,
            tgtPatch,
            tgtToSrc
        );

    // Make sure every face references exactly one face
    auto forwardCheck = []
    (
        const primitiveOldTimePatch& patch,
        const List<DynamicList<label>>& localOtherFaces,
        const bool isSrc
    )
    {
        forAll(localOtherFaces, facei)
        {
            if (localOtherFaces[facei].size() != 1)
            {
                FatalErrorInFunction
                    << (isSrc ? "Source" : "Target")
                    << " face #" << facei << " at "
                    << patch.faceCentres()[facei]
                    << " did not match a face on the "
                    << (isSrc ? "target" : "source")
                    << " side" << exit(FatalError);
            }
        }
    };
    forwardCheck(srcPatch, srcLocalTgtFaces_, true);
    forwardCheck(tgtPatch, tgtLocalSrcFaces_, false);

    // Make sure every face is referenced by exactly one face
    auto reverseCheck = []
    (
        const primitiveOldTimePatch& patch,
        const List<DynamicList<label>>& otherLocalFaces,
        const autoPtr<distributionMap>& mapPtr,
        const bool isSrc
    )
    {
        labelList count
        (
            mapPtr.valid() ? mapPtr->constructSize() : patch.size(),
            0
        );

        forAll(otherLocalFaces, otherFacei)
        {
            forAll(otherLocalFaces[otherFacei], i)
            {
                count[otherLocalFaces[otherFacei][i]] ++;
            }
        }

        if (mapPtr.valid())
        {
            distributionMapBase::distribute
            (
                Pstream::commsTypes::nonBlocking,
                List<labelPair>(),
                patch.size(),
                mapPtr->constructMap(),
                false,
                mapPtr->subMap(),
                false,
                count,
                plusEqOp<label>(),
                flipOp(),
                label(0)
            );
        }

        forAll(count, facei)
        {
            if (count[facei] != 1)
            {
                FatalErrorInFunction
                    << (isSrc ? "Source" : "Target")
                    << " face #" << facei << " at "
                    << patch.faceCentres()[facei]
                    << " did not match a face on the "
                    << (isSrc ? "target" : "source")
                    << " side" << exit(FatalError);
            }
        }
    };
    reverseCheck(srcPatch, tgtLocalSrcFaces_, srcMapPtr_, true);
    reverseCheck(tgtPatch, srcLocalTgtFaces_, tgtMapPtr_, false);

    return nCouples;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatches::matching::matching(const bool reverse)
:
    nearest(reverse)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToPatches::matching::~matching()
{}


// ************************************************************************* //
