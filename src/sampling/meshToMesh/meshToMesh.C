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
#include "processorPolyPatch.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshToMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::meshToMesh::calculate(const word& methodName)
{
    Info<< "Creating mesh-to-mesh addressing for " << srcRegion_.name()
        << " and " << tgtRegion_.name() << " regions using "
        << methodName << endl;

    singleMeshProc_ = calcDistribution(srcRegion_, tgtRegion_);

    autoPtr<meshToMeshMethod> methodPtr;

    scalar V = 0;

    if (singleMeshProc_ == -1)
    {
        // create global indexing for src and tgt meshes
        globalIndex globalSrcCells(srcRegion_.nCells());
        globalIndex globalTgtCells(tgtRegion_.nCells());

        // Create processor map of overlapping cells. This map gets
        // (possibly remote) cells from the tgt mesh such that they (together)
        // cover all of the src mesh
        autoPtr<distributionMap> mapPtr = calcProcMap(srcRegion_, tgtRegion_);
        const distributionMap& map = mapPtr();

        // Distribute the target geometry to create local-target geometry
        // (i.e., local to the source)
        pointField localTgtPoints;
        faceList localTgtFaces;
        labelList localTgtFaceOwners;
        labelList localTgtFaceNeighbours;
        labelList localTgtGlobalCellIDs;
        distributeAndMergeCells
        (
            map,
            tgtRegion_,
            globalTgtCells,
            localTgtPoints,
            localTgtFaces,
            localTgtFaceOwners,
            localTgtFaceNeighbours,
            localTgtGlobalCellIDs
        );

        // create a new target mesh
        polyMesh localTgtMesh
        (
            IOobject
            (
                "localTgt",
                tgtRegion_.time().timeName(),
                tgtRegion_.time(),
                IOobject::NO_READ
            ),
            move(localTgtPoints),
            move(localTgtFaces),
            move(localTgtFaceOwners),
            move(localTgtFaceNeighbours),
            false
        );

        // create some dummy patch info and add to the new target mesh
        List<polyPatch*> patches(1);
        patches[0] = new polyPatch
        (
            "defaultFaces",
            localTgtMesh.nFaces() - localTgtMesh.nInternalFaces(),
            localTgtMesh.nInternalFaces(),
            0,
            localTgtMesh.boundaryMesh(),
            word::null
        );
        localTgtMesh.addPatches(patches);

        // force calculation of tet-base points used for point-in-cell
        (void)localTgtMesh.tetBasePtIs();

        if (debug)
        {
            Pout<< "Created local target mesh:" << nl
                << " cells = " << tgtRegion_.nCells()
                << ", local cells = " << localTgtMesh.nCells() << nl
                << " faces = " << tgtRegion_.nFaces()
                << ", local faces = " << localTgtMesh.nFaces() << endl;

            if (debug > 1)
            {
                Pout<< "Writing local target mesh: "
                    << localTgtMesh.name() << endl;
                localTgtMesh.write();
            }
        }

        // Do the intersection
        methodPtr = meshToMeshMethod::New(methodName, srcRegion_, localTgtMesh);
        V = methodPtr->calculate
            (
                srcToTgtCellAddr_,
                srcToTgtCellWght_,
                tgtToSrcCellAddr_,
                tgtToSrcCellWght_
            );

        // per source cell the target cell address in localTgtMesh mesh
        forAll(srcToTgtCellAddr_, i)
        {
            labelList& addressing = srcToTgtCellAddr_[i];
            forAll(addressing, addrI)
            {
                addressing[addrI] = localTgtGlobalCellIDs[addressing[addrI]];
            }
        }

        // convert target addresses in localTgtMesh into global cell numbering
        forAll(tgtToSrcCellAddr_, i)
        {
            labelList& addressing = tgtToSrcCellAddr_[i];
            forAll(addressing, addrI)
            {
                addressing[addrI] = globalSrcCells.toGlobal(addressing[addrI]);
            }
        }

        // send the connectivity back to the target
        distributionMapBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            tgtRegion_.nCells(),
            map.constructMap(),
            false,
            map.subMap(),
            false,
            tgtToSrcCellAddr_,
            ListAppendEqOp<label>(),
            flipOp(),
            labelList()
        );
        distributionMapBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            tgtRegion_.nCells(),
            map.constructMap(),
            false,
            map.subMap(),
            false,
            tgtToSrcCellWght_,
            ListAppendEqOp<scalar>(),
            flipOp(),
            scalarList()
        );

        // weights normalisation
        methodPtr->normalise
        (
            srcRegion_,
            srcToTgtCellAddr_,
            srcToTgtCellWght_
        );
        methodPtr->normalise
        (
            tgtRegion_,
            tgtToSrcCellAddr_,
            tgtToSrcCellWght_
        );

        // cache maps and reset addresses
        List<Map<label>> cMap;
        srcMapPtr_.reset
        (
            new distributionMap(globalSrcCells, tgtToSrcCellAddr_, cMap)
        );
        tgtMapPtr_.reset
        (
            new distributionMap(globalTgtCells, srcToTgtCellAddr_, cMap)
        );

        // collect volume intersection contributions
        reduce(V, sumOp<scalar>());
    }
    else
    {
        // Do the intersection
        methodPtr = meshToMeshMethod::New(methodName, srcRegion_, tgtRegion_);
        V = methodPtr->calculate
            (
                srcToTgtCellAddr_,
                srcToTgtCellWght_,
                tgtToSrcCellAddr_,
                tgtToSrcCellWght_
            );

        // Normalise the weights
        methodPtr->normalise
        (
            srcRegion_,
            srcToTgtCellAddr_,
            srcToTgtCellWght_
        );
        methodPtr->normalise
        (
            tgtRegion_,
            tgtToSrcCellAddr_,
            tgtToSrcCellWght_
        );
    }

    Info<< "    Overlap volume: " << V << endl;

    return V;
}


void Foam::meshToMesh::calculatePatchToPatches(const word& methodName)
{
    if (!patchToPatches_.empty())
    {
        FatalErrorInFunction
            << "patchToPatches already calculated"
            << exit(FatalError);
    }

    const word& patchToPatchType =
        meshToMeshMethod::New
        (
            methodName,
            NullObjectRef<polyMesh>(),
            NullObjectRef<polyMesh>()
        )->patchToPatchMethod();

    patchToPatches_.setSize(srcPatchID_.size());

    forAll(srcPatchID_, i)
    {
        const label srcPatchi = srcPatchID_[i];
        const label tgtPatchi = tgtPatchID_[i];

        const polyPatch& srcPP = srcRegion_.boundaryMesh()[srcPatchi];
        const polyPatch& tgtPP = tgtRegion_.boundaryMesh()[tgtPatchi];

        Info<< "Creating patchToPatch between source patch "
            << srcPP.name() << " and target patch " << tgtPP.name()
            << " using " << patchToPatchType << endl;

        Info<< incrIndent;

        patchToPatches_.set(i, patchToPatch::New(patchToPatchType, true));

        patchToPatches_[i].update
        (
            srcPP,
            PatchTools::pointNormals(srcRegion_, srcPP),
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
        const polyBoundaryMesh& srcBM = srcRegion_.boundaryMesh();
        const polyBoundaryMesh& tgtBM = tgtRegion_.boundaryMesh();

        DynamicList<label> srcPatchID(srcBM.size());
        DynamicList<label> tgtPatchID(tgtBM.size());
        forAll(srcBM, patchi)
        {
            const polyPatch& pp = srcBM[patchi];

            // We want to map all the global patches, including constraint
            // patches (since they might have mappable properties, e.g.
            // jumpCyclic). We'll fix the value afterwards.
            if (!isA<processorPolyPatch>(pp))
            {
                srcPatchID.append(pp.index());

                label tgtPatchi = tgtBM.findPatchID(pp.name());

                if (tgtPatchi != -1)
                {
                    tgtPatchID.append(tgtPatchi);
                }
                else
                {
                    FatalErrorInFunction
                        << "Source patch " << pp.name()
                        << " not found in target mesh. "
                        << "Available target patches are " << tgtBM.names()
                        << exit(FatalError);
                }
            }
        }

        srcPatchID_.transfer(srcPatchID);
        tgtPatchID_.transfer(tgtPatchID);
    }

    // calculate volume addressing and weights
    calculate(methodName);

    // calculate patch addressing and weights
    calculatePatchToPatches(methodName);
}


void Foam::meshToMesh::constructFromCuttingPatches
(
    const word& methodName,
    const HashTable<word>& patchMap,
    const wordList& cuttingPatches
)
{
    srcPatchID_.setSize(patchMap.size());
    tgtPatchID_.setSize(patchMap.size());

    label i = 0;
    forAllConstIter(HashTable<word>, patchMap, iter)
    {
        const word& tgtPatchName = iter.key();
        const word& srcPatchName = iter();

        const polyPatch& srcPatch = srcRegion_.boundaryMesh()[srcPatchName];
        const polyPatch& tgtPatch = tgtRegion_.boundaryMesh()[tgtPatchName];

        srcPatchID_[i] = srcPatch.index();
        tgtPatchID_[i] = tgtPatch.index();
        i++;
    }

    // calculate volume addressing and weights
    calculate(methodName);

    // calculate patch addressing and weights
    calculatePatchToPatches(methodName);

    // set IDs of cutting patches on target mesh
    cuttingPatches_.setSize(cuttingPatches.size());
    forAll(cuttingPatches_, i)
    {
        const word& patchName = cuttingPatches[i];
        cuttingPatches_[i] = tgtRegion_.boundaryMesh().findPatchID(patchName);
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
    srcRegion_(src),
    tgtRegion_(tgt),
    srcPatchID_(),
    tgtPatchID_(),
    patchToPatches_(),
    cuttingPatches_(),
    srcToTgtCellAddr_(),
    tgtToSrcCellAddr_(),
    srcToTgtCellWght_(),
    tgtToSrcCellWght_(),
    singleMeshProc_(-1),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr)
{
    constructNoCuttingPatches(methodName, interpAllPatches);
}


Foam::meshToMesh::meshToMesh
(
    const polyMesh& src,
    const polyMesh& tgt,
    const word& methodName,
    const HashTable<word>& patchMap,
    const wordList& cuttingPatches
)
:
    srcRegion_(src),
    tgtRegion_(tgt),
    srcPatchID_(),
    tgtPatchID_(),
    patchToPatches_(),
    cuttingPatches_(),
    srcToTgtCellAddr_(),
    tgtToSrcCellAddr_(),
    srcToTgtCellWght_(),
    tgtToSrcCellWght_(),
    singleMeshProc_(-1),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr)
{
    constructFromCuttingPatches
    (
        methodName,
        patchMap,
        cuttingPatches
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshToMesh::~meshToMesh()
{}


// ************************************************************************* //
