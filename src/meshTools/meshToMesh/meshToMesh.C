/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2026 OpenFOAM Foundation
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

#include "meshToMesh.H"
#include "PatchTools.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "processorPolyPatch.H"
#include "nearestCellsToCellsExtrapolation.H"
#include "nearestPatchToPatchExtrapolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshToMesh, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshToMesh::meshToMesh
(
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh,
    const word& cellsInterpolationType,
    const word& cellsExtrapolationType,
    const word& patchInterpolationType,
    const word& patchExtrapolationType,
    const HashTable<word>& patchMap
)
:
    srcMesh_(srcMesh),
    tgtMesh_(tgtMesh),
    cellsInterpolation_(),
    cellsExtrapolationType_(cellsExtrapolationType),
    srcCellsExtrapolation_(),
    tgtCellsExtrapolation_(),
    patchIndices_(),
    patchInterpolations_(),
    patchExtrapolationType_(patchExtrapolationType),
    srcPatchExtrapolations_(),
    tgtPatchExtrapolations_()
{
    // If no patch map was supplied, then assume a consistent pair of meshes in
    // which corresponding patches have the same name
    if (isNull(patchMap))
    {
        DynamicList<labelPair> patchIDs;

        forAll(srcMesh_.boundary(), srcPatchi)
        {
            const polyPatch& srcPp = srcMesh_.boundary()[srcPatchi];

            // We don't want to map empty or wedge patches, as by definition
            // these do not hold relevant values. We also don't want to map
            // processor patches as these are likely to differ between cases.
            // In general, though, we do want to map constraint patches as they
            // might have additional mappable properties; e.g., the jump field
            // of a jump cyclic.
            if
            (
                !isA<emptyPolyPatch>(srcPp)
             && !isA<wedgePolyPatch>(srcPp)
             && !isA<processorPolyPatch>(srcPp)
            )
            {
                const label tgtPatchi =
                    tgtMesh_.boundary().findIndex(srcPp.name());

                if (tgtPatchi == -1)
                {
                    FatalErrorInFunction
                        << "Source patch " << srcPp.name()
                        << " not found in target mesh. "
                        << "Available target patches are "
                        << tgtMesh_.boundary().names()
                        << exit(FatalError);
                }

                patchIDs.append(labelPair(srcPatchi, tgtPatchi));
            }
        }

        patchIndices_.transfer(patchIDs);
    }

    // If a patch mas was supplied then convert it to pairs of patch indices
    else
    {
        patchIndices_.setSize(patchMap.size());
        label i = 0;
        forAllConstIter(HashTable<word>, patchMap, iter)
        {
            const word& tgtPatchName = iter.key();
            const word& srcPatchName = iter();

            const label srcPatchi =
                srcMesh_.boundary().findIndex(srcPatchName);
            const label tgtPatchi =
                tgtMesh_.boundary().findIndex(tgtPatchName);

            if (srcPatchi == -1)
            {
                FatalErrorInFunction
                    << "Patch " << srcPatchName
                    << " not found in source mesh. "
                    << "Available source patches are "
                    << srcMesh_.boundary().names()
                    << exit(FatalError);
            }
            if (tgtPatchi == -1)
            {
                FatalErrorInFunction
                    << "Patch " << tgtPatchName
                    << " not found in target mesh. "
                    << "Available target patches are "
                    << tgtMesh_.boundary().names()
                    << exit(FatalError);
            }

            patchIndices_[i ++] = labelPair(srcPatchi, tgtPatchi);
        }
    }

    // Calculate cell addressing and weights
    Info<< "Creating cellsToCells between source mesh "
        << srcMesh_.name() << " and target mesh " << tgtMesh_.name()
        << " using " << cellsInterpolationType << endl << incrIndent;

    cellsInterpolation_ = cellsToCells::New(cellsInterpolationType);
    cellsInterpolation_->update(srcMesh_, tgtMesh_);

    srcCellsExtrapolation_.clear();
    tgtCellsExtrapolation_.clear();

    Info<< decrIndent;

    // Calculate patch addressing and weights
    patchInterpolations_.setSize(patchIndices_.size());
    srcPatchExtrapolations_.setSize(patchIndices_.size());
    tgtPatchExtrapolations_.setSize(patchIndices_.size());
    forAll(patchIndices_, i)
    {
        const label srcPatchi = patchIndices_[i].first();
        const label tgtPatchi = patchIndices_[i].second();

        const polyPatch& srcPp = srcMesh_.boundary()[srcPatchi];
        const polyPatch& tgtPp = tgtMesh_.boundary()[tgtPatchi];

        Info<< "Creating patchToPatch between source patch "
            << srcPp.name() << " and target patch " << tgtPp.name()
            << " using " << patchInterpolationType << endl << incrIndent;

        patchInterpolations_.set
        (
            i,
            patchToPatch::New(patchInterpolationType, true)
        );

        patchInterpolations_[i].update
        (
            srcPp,
            PatchTools::pointNormals(srcMesh_, srcPp),
            tgtPp
        );

        Info<< decrIndent;
    }
}


Foam::meshToMesh::meshToMesh
(
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh,
    const word& cellsInterpolationType,
    const word& patchInterpolationType,
    const HashTable<word>& patchMap
)
:
    meshToMesh
    (
        srcMesh,
        tgtMesh,
        cellsInterpolationType,
        cellsToCellsExtrapolations::nearest::typeName,
        patchInterpolationType,
        patchToPatchExtrapolations::nearest::typeName,
        patchMap
    )
{}


Foam::meshToMesh::meshToMesh
(
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh,
    const word& interpolationType,
    const HashTable<word>& patchMap
)
:
    meshToMesh
    (
        srcMesh,
        tgtMesh,
        interpolationType,
        cellsToCellsExtrapolations::nearest::typeName,
        interpolationType,
        patchToPatchExtrapolations::nearest::typeName,
        patchMap
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshToMesh::~meshToMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::meshToMesh::consistent() const
{
    boolList srcPatchIsMapped(srcMesh_.boundary().size(), false);
    boolList tgtPatchIsMapped(tgtMesh_.boundary().size(), false);

    // Mark anything paired as mapped
    forAll(patchIndices_, i)
    {
        const label srcPatchi = patchIndices_[i].first();
        const label tgtPatchi = patchIndices_[i].second();

        srcPatchIsMapped[srcPatchi] = true;
        tgtPatchIsMapped[tgtPatchi] = true;
    }

    // Filter out un-mappable patches
    forAll(srcMesh_.boundary(), srcPatchi)
    {
        const polyPatch& srcPp = srcMesh_.boundary()[srcPatchi];

        if
        (
            isA<emptyPolyPatch>(srcPp)
         || isA<wedgePolyPatch>(srcPp)
         || isA<processorPolyPatch>(srcPp)
        )
        {
            srcPatchIsMapped[srcPp.index()] = true;
        }
    }
    forAll(tgtMesh_.boundary(), tgtPatchi)
    {
        const polyPatch& tgtPp = tgtMesh_.boundary()[tgtPatchi];

        if
        (
            isA<emptyPolyPatch>(tgtPp)
         || isA<wedgePolyPatch>(tgtPp)
         || isA<processorPolyPatch>(tgtPp)
        )
        {
            tgtPatchIsMapped[tgtPp.index()] = true;
        }
    }

    // Return whether or not everything is mapped
    return
        findIndex(srcPatchIsMapped, false) == -1
     && findIndex(tgtPatchIsMapped, false) == -1;
}


Foam::remote Foam::meshToMesh::srcToTgtPoint
(
    const label srcCelli,
    const point& p
) const
{
    return cellsInterpolation_().srcToTgtPoint(tgtMesh_, srcCelli, p);
}


// ************************************************************************* //
