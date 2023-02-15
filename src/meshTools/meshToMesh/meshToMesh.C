/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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
    const word& engineType,
    const HashTable<word>& patchMap
)
:
    srcMesh_(srcMesh),
    tgtMesh_(tgtMesh),
    cellsInterpolation_(),
    srcCellsStabilisation_(),
    tgtCellsStabilisation_(),
    patchIDs_(),
    patchInterpolations_(),
    srcPatchStabilisations_(),
    tgtPatchStabilisations_()
{
    // If no patch map was supplied, then assume a consistent pair of meshes in
    // which corresponding patches have the same name
    if (isNull(patchMap))
    {
        DynamicList<labelPair> patchIDs;

        forAll(srcMesh_.boundaryMesh(), srcPatchi)
        {
            const polyPatch& srcPp = srcMesh_.boundaryMesh()[srcPatchi];

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
                    tgtMesh_.boundaryMesh().findPatchID(srcPp.name());

                if (tgtPatchi == -1)
                {
                    FatalErrorInFunction
                        << "Source patch " << srcPp.name()
                        << " not found in target mesh. "
                        << "Available target patches are "
                        << tgtMesh_.boundaryMesh().names()
                        << exit(FatalError);
                }

                patchIDs.append(labelPair(srcPatchi, tgtPatchi));
            }
        }

        patchIDs_.transfer(patchIDs);
    }

    // If a patch mas was supplied then convert it to pairs of patch indices
    else
    {
        patchIDs_.setSize(patchMap.size());
        label i = 0;
        forAllConstIter(HashTable<word>, patchMap, iter)
        {
            const word& tgtPatchName = iter.key();
            const word& srcPatchName = iter();

            const label srcPatchi =
                srcMesh_.boundaryMesh().findPatchID(srcPatchName);
            const label tgtPatchi =
                tgtMesh_.boundaryMesh().findPatchID(tgtPatchName);

            if (srcPatchi == -1)
            {
                FatalErrorInFunction
                    << "Patch " << srcPatchName
                    << " not found in source mesh. "
                    << "Available source patches are "
                    << srcMesh_.boundaryMesh().names()
                    << exit(FatalError);
            }
            if (tgtPatchi == -1)
            {
                FatalErrorInFunction
                    << "Patch " << tgtPatchName
                    << " not found in target mesh. "
                    << "Available target patches are "
                    << tgtMesh_.boundaryMesh().names()
                    << exit(FatalError);
            }

            patchIDs_[i ++] = labelPair(srcPatchi, tgtPatchi);
        }
    }

    // Calculate cell addressing and weights
    Info<< "Creating cellsToCells between source mesh "
        << srcMesh_.name() << " and target mesh " << tgtMesh_.name()
        << " using " << engineType << endl << incrIndent;

    cellsInterpolation_ = cellsToCells::New(engineType);
    cellsInterpolation_->update(srcMesh_, tgtMesh_);

    srcCellsStabilisation_.clear();
    tgtCellsStabilisation_.clear();

    Info<< decrIndent;

    // Calculate patch addressing and weights
    patchInterpolations_.setSize(patchIDs_.size());
    srcPatchStabilisations_.setSize(patchIDs_.size());
    tgtPatchStabilisations_.setSize(patchIDs_.size());
    forAll(patchIDs_, i)
    {
        const label srcPatchi = patchIDs_[i].first();
        const label tgtPatchi = patchIDs_[i].second();

        const polyPatch& srcPp = srcMesh_.boundaryMesh()[srcPatchi];
        const polyPatch& tgtPp = tgtMesh_.boundaryMesh()[tgtPatchi];

        Info<< "Creating patchToPatch between source patch "
            << srcPp.name() << " and target patch " << tgtPp.name()
            << " using " << engineType << endl << incrIndent;

        patchInterpolations_.set
        (
            i,
            patchToPatch::New(engineType, true)
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


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshToMesh::~meshToMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::meshToMesh::consistent() const
{
    boolList srcPatchIsMapped(srcMesh_.boundaryMesh().size(), false);
    boolList tgtPatchIsMapped(tgtMesh_.boundaryMesh().size(), false);

    // Mark anything paired as mapped
    forAll(patchIDs_, i)
    {
        const label srcPatchi = patchIDs_[i].first();
        const label tgtPatchi = patchIDs_[i].second();

        srcPatchIsMapped[srcPatchi] = true;
        tgtPatchIsMapped[tgtPatchi] = true;
    }

    // Filter out un-mappable patches
    forAll(srcMesh_.boundaryMesh(), srcPatchi)
    {
        const polyPatch& srcPp = srcMesh_.boundaryMesh()[srcPatchi];

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
    forAll(tgtMesh_.boundaryMesh(), tgtPatchi)
    {
        const polyPatch& tgtPp = tgtMesh_.boundaryMesh()[tgtPatchi];

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
