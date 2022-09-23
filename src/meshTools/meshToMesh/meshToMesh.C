/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2022 OpenFOAM Foundation
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
#include "globalIndex.H"
#include "meshToMeshMethod.H"
#include "PatchTools.H"
#include "patchToPatchTools.H"
#include "processorPolyPatch.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshToMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::meshToMesh::calculateCellToCells(const word& methodName)
{
    Info<< "Creating mesh-to-mesh addressing for " << srcMesh_.name()
        << " and " << tgtMesh_.name() << " regions using "
        << methodName << endl;

    singleProcess_ =
        patchToPatchTools::singleProcess
        (
            srcMesh_.nCells(),
            tgtMesh_.nCells()
        );

    autoPtr<meshToMeshMethod> methodPtr;

    scalar V = 0;

    if (isSingleProcess())
    {
        // Do the intersection
        methodPtr = meshToMeshMethod::New(methodName, srcMesh_, tgtMesh_);
        V = methodPtr->calculate
            (
                srcLocalTgtCells_,
                srcWeights_,
                tgtLocalSrcCells_,
                tgtWeights_
            );

        // Normalise the weights
        methodPtr->normalise(srcMesh_, srcLocalTgtCells_, srcWeights_);
        methodPtr->normalise(tgtMesh_, tgtLocalSrcCells_, tgtWeights_);
    }
    else
    {
        // Create the target map of overlapping cells. This map gets remote
        // parts of the target mesh so that everything needed to compute an
        // intersection is available locally to the source. Use it to create a
        // source-local target mesh.
        tgtMapPtr_ =
            patchToPatchTools::constructDistributionMap
            (
                tgtMeshSendCells(srcMesh_, tgtMesh_)
            );
        localTgtProcCellsPtr_.reset
        (
            new List<remote>
            (
                distributeMesh
                (
                    tgtMapPtr_(),
                    tgtMesh_,
                    localTgtMeshPtr_
                )
            )
        );
        const polyMesh& localTgtMesh = localTgtMeshPtr_();

        if (debug > 1)
        {
            Pout<< "Writing local target mesh: "
                << localTgtMesh.name() << endl;
            localTgtMesh.write();
        }

        // Do the intersection
        methodPtr = meshToMeshMethod::New(methodName, srcMesh_, localTgtMesh);
        V = methodPtr->calculate
            (
                srcLocalTgtCells_,
                srcWeights_,
                tgtLocalSrcCells_,
                tgtWeights_
            );

        // Trim the local target mesh
        trimLocalTgt();

        if (debug > 1)
        {
            Pout<< "Writing trimmed local target mesh: "
                << localTgtMesh.name() << endl;
            localTgtMesh.write();
        }

        // Construct the source map
        srcMapPtr_ =
            patchToPatchTools::constructDistributionMap
            (
                patchToPatchTools::procSendIndices
                (
                    tgtLocalSrcCells_,
                    localTgtProcCellsPtr_()
                )
            );
        localSrcProcCellsPtr_.reset
        (
            new List<remote>
            (
                patchToPatchTools::distributeAddressing(srcMapPtr_())
            )
        );

        // Collect the addressing on the target
        patchToPatchTools::rDistributeTgtAddressing
        (
            tgtMesh_.nCells(),
            tgtMapPtr_(),
            localSrcProcCellsPtr_(),
            tgtLocalSrcCells_
        );

        // Collect the weights on the target
        patchToPatchTools::rDistributeListList
        (
            tgtMesh_.nCells(),
            tgtMapPtr_(),
            tgtWeights_
        );

        // Normalise the weights
        methodPtr->normalise(srcMesh_, srcLocalTgtCells_, srcWeights_);
        methodPtr->normalise(tgtMesh_, tgtLocalSrcCells_, tgtWeights_);

        // collect volume intersection contributions
        reduce(V, sumOp<scalar>());
    }

    Info<< "    Overlap volume: " << V << endl;

    return V;
}


void Foam::meshToMesh::calculatePatchToPatches(const word& methodName)
{
    if (!srcToTgtPatchToPatches_.empty())
    {
        FatalErrorInFunction
            << "srcToTgtPatchToPatches already calculated"
            << exit(FatalError);
    }

    const word& patchToPatchType =
        meshToMeshMethod::New
        (
            methodName,
            NullObjectRef<polyMesh>(),
            NullObjectRef<polyMesh>()
        )->patchToPatchMethod();

    srcToTgtPatchToPatches_.setSize(srcToTgtPatchIDs_.size());

    forAll(srcToTgtPatchIDs_, i)
    {
        const label srcPatchi = srcToTgtPatchIDs_[i].first();
        const label tgtPatchi = srcToTgtPatchIDs_[i].second();

        const polyPatch& srcPP = srcMesh_.boundaryMesh()[srcPatchi];
        const polyPatch& tgtPP = tgtMesh_.boundaryMesh()[tgtPatchi];

        Info<< "Creating patchToPatch between source patch "
            << srcPP.name() << " and target patch " << tgtPP.name()
            << " using " << patchToPatchType << endl;

        Info<< incrIndent;

        srcToTgtPatchToPatches_.set
        (
            i,
            patchToPatch::New(patchToPatchType, true)
        );

        srcToTgtPatchToPatches_[i].update
        (
            srcPP,
            PatchTools::pointNormals(srcMesh_, srcPP),
            tgtPP
        );

        Info<< decrIndent;
    }
}


void Foam::meshToMesh::constructNoCuttingPatches
(
    const word& methodName,
    const bool interpAllPatches
)
{
    if (interpAllPatches)
    {
        DynamicList<labelPair> srcToTgtPatchIDs;

        forAll(srcMesh_.boundaryMesh(), srcPatchi)
        {
            const polyPatch& srcPp = srcMesh_.boundaryMesh()[srcPatchi];

            // We want to map all the global patches, including constraint
            // patches (since they might have mappable properties, e.g.
            // jumpCyclic). We'll fix the value afterwards.
            if (!isA<processorPolyPatch>(srcPp))
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

                srcToTgtPatchIDs.append(labelPair(srcPatchi, tgtPatchi));
            }
        }

        srcToTgtPatchIDs_.transfer(srcToTgtPatchIDs);
    }

    // calculate patch addressing and weights
    calculatePatchToPatches(methodName);

    // calculate cell addressing and weights
    calculateCellToCells(methodName);
}


void Foam::meshToMesh::constructFromCuttingPatches
(
    const word& methodName,
    const HashTable<word>& patchMap,
    const wordList& tgtCuttingPatches
)
{
    srcToTgtPatchIDs_.setSize(patchMap.size());
    label i = 0;
    forAllConstIter(HashTable<word>, patchMap, iter)
    {
        const word& tgtPatchName = iter.key();
        const word& srcPatchName = iter();

        srcToTgtPatchIDs_[i++] =
            labelPair
            (
                srcMesh_.boundaryMesh().findPatchID(srcPatchName),
                tgtMesh_.boundaryMesh().findPatchID(tgtPatchName)
            );
    }

    // calculate patch addressing and weights
    calculatePatchToPatches(methodName);

    // calculate cell addressing and weights
    calculateCellToCells(methodName);

    // set IDs of cutting patches on target mesh
    tgtCuttingPatchIDs_.setSize(tgtCuttingPatches.size());
    forAll(tgtCuttingPatchIDs_, i)
    {
        const word& patchName = tgtCuttingPatches[i];
        tgtCuttingPatchIDs_[i] = tgtMesh_.boundaryMesh().findPatchID(patchName);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshToMesh::meshToMesh
(
    const polyMesh& src,
    const polyMesh& tgt,
    const word& methodName,
    bool interpAllPatches
)
:
    srcMesh_(src),
    tgtMesh_(tgt),
    srcToTgtPatchIDs_(),
    srcToTgtPatchToPatches_(),
    tgtCuttingPatchIDs_(),
    srcLocalTgtCells_(),
    tgtLocalSrcCells_(),
    srcWeights_(),
    tgtWeights_(),
    singleProcess_(-1),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr),
    localSrcProcCellsPtr_(nullptr),
    localTgtProcCellsPtr_(nullptr),
    localTgtMeshPtr_(nullptr)
{
    constructNoCuttingPatches
    (
        methodName,
        interpAllPatches
    );
}


Foam::meshToMesh::meshToMesh
(
    const polyMesh& src,
    const polyMesh& tgt,
    const word& methodName,
    const HashTable<word>& patchMap,
    const wordList& tgtCuttingPatches
)
:
    srcMesh_(src),
    tgtMesh_(tgt),
    srcToTgtPatchIDs_(),
    srcToTgtPatchToPatches_(),
    tgtCuttingPatchIDs_(),
    srcLocalTgtCells_(),
    tgtLocalSrcCells_(),
    srcWeights_(),
    tgtWeights_(),
    singleProcess_(-1),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr),
    localTgtMeshPtr_(nullptr)
{
    constructFromCuttingPatches
    (
        methodName,
        patchMap,
        tgtCuttingPatches
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshToMesh::~meshToMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::remote Foam::meshToMesh::srcToTgtPoint
(
    const label srcCelli,
    const point& p
) const
{
    forAll(srcLocalTgtCells_[srcCelli], i)
    {
        const label tgtCelli = srcLocalTgtCells_[srcCelli][i];

        const polyMesh& tgtMesh =
            singleProcess_ == -1 ? localTgtMeshPtr_() : tgtMesh_;

        if (tgtMesh.pointInCell(p, tgtCelli))
        {
            return
                singleProcess_ == -1
              ? localTgtProcCellsPtr_()[tgtCelli]
              : remote(Pstream::myProcNo(), tgtCelli);
        }
    }

    return remote();
}


// ************************************************************************* //
