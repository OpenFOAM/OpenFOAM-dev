/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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
#include "Time.H"
#include "globalIndex.H"
#include "meshToMeshMethod.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshToMesh, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::meshToMesh::interpolationMethod,
        3
    >::names[] =
    {
        "direct",
        "mapNearest",
        "cellVolumeWeight"
    };

    const NamedEnum<meshToMesh::interpolationMethod, 3>
        meshToMesh::interpolationMethodNames_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<>
void Foam::meshToMesh::mapAndOpSrcToTgt
(
    const AMIInterpolation& AMI,
    const Field<scalar>& srcField,
    Field<scalar>& tgtField,
    const plusEqOp<scalar>& cop
) const
{}


template<>
void Foam::meshToMesh::mapAndOpSrcToTgt
(
    const AMIInterpolation& AMI,
    const Field<vector>& srcField,
    Field<vector>& tgtField,
    const plusEqOp<vector>& cop
) const
{}


template<>
void Foam::meshToMesh::mapAndOpSrcToTgt
(
    const AMIInterpolation& AMI,
    const Field<sphericalTensor>& srcField,
    Field<sphericalTensor>& tgtField,
    const plusEqOp<sphericalTensor>& cop
) const
{}


template<>
void Foam::meshToMesh::mapAndOpSrcToTgt
(
    const AMIInterpolation& AMI,
    const Field<symmTensor>& srcField,
    Field<symmTensor>& tgtField,
    const plusEqOp<symmTensor>& cop
) const
{}


template<>
void Foam::meshToMesh::mapAndOpSrcToTgt
(
    const AMIInterpolation& AMI,
    const Field<tensor>& srcField,
    Field<tensor>& tgtField,
    const plusEqOp<tensor>& cop
) const
{}


template<>
void Foam::meshToMesh::mapAndOpTgtToSrc
(
    const AMIInterpolation& AMI,
    Field<scalar>& srcField,
    const Field<scalar>& tgtField,
    const plusEqOp<scalar>& cop
) const
{}


template<>
void Foam::meshToMesh::mapAndOpTgtToSrc
(
    const AMIInterpolation& AMI,
    Field<vector>& srcField,
    const Field<vector>& tgtField,
    const plusEqOp<vector>& cop
) const
{}


template<>
void Foam::meshToMesh::mapAndOpTgtToSrc
(
    const AMIInterpolation& AMI,
    Field<sphericalTensor>& srcField,
    const Field<sphericalTensor>& tgtField,
    const plusEqOp<sphericalTensor>& cop
) const
{}


template<>
void Foam::meshToMesh::mapAndOpTgtToSrc
(
    const AMIInterpolation& AMI,
    Field<symmTensor>& srcField,
    const Field<symmTensor>& tgtField,
    const plusEqOp<symmTensor>& cop
) const
{}


template<>
void Foam::meshToMesh::mapAndOpTgtToSrc
(
    const AMIInterpolation& AMI,
    Field<tensor>& srcField,
    const Field<tensor>& tgtField,
    const plusEqOp<tensor>& cop
) const
{}


Foam::labelList Foam::meshToMesh::maskCells
(
    const polyMesh& src,
    const polyMesh& tgt
) const
{
    boundBox intersectBb
    (
        max(src.bounds().min(), tgt.bounds().min()),
        min(src.bounds().max(), tgt.bounds().max())
    );

    intersectBb.inflate(0.01);

    const cellList& srcCells = src.cells();
    const faceList& srcFaces = src.faces();
    const pointField& srcPts = src.points();

    DynamicList<label> cells(src.size());
    forAll(srcCells, srcI)
    {
        boundBox cellBb(srcCells[srcI].points(srcFaces, srcPts), false);
        if (intersectBb.overlaps(cellBb))
        {
            cells.append(srcI);
        }
    }

    if (debug)
    {
        Pout<< "participating source mesh cells: " << cells.size() << endl;
    }

    return move(cells);
}


void Foam::meshToMesh::normaliseWeights
(
    const word& descriptor,
    const labelListList& addr,
    scalarListList& wght
) const
{
    const label nCell = returnReduce(wght.size(), sumOp<label>());

    if (nCell > 0)
    {
        forAll(wght, celli)
        {
            scalarList& w = wght[celli];
            scalar s = sum(w);

            forAll(w, i)
            {
                // note: normalise by s instead of cell volume since
                // 1-to-1 methods duplicate contributions in parallel
                w[i] /= s;
            }
        }
    }
}


void Foam::meshToMesh::calcAddressing
(
    const word& methodName,
    const polyMesh& src,
    const polyMesh& tgt
)
{
    autoPtr<meshToMeshMethod> methodPtr
    (
        meshToMeshMethod::New
        (
            methodName,
            src,
            tgt
        )
    );

    methodPtr->calculate
    (
        srcToTgtCellAddr_,
        srcToTgtCellWght_,
        tgtToSrcCellAddr_,
        tgtToSrcCellWght_
    );

    V_ = methodPtr->V();

    if (debug)
    {
        methodPtr->writeConnectivity(src, tgt, srcToTgtCellAddr_);
    }
}


void Foam::meshToMesh::calculate(const word& methodName)
{
    Info<< "Creating mesh-to-mesh addressing for " << srcRegion_.name()
        << " and " << tgtRegion_.name() << " regions using "
        << methodName << endl;

    singleMeshProc_ = calcDistribution(srcRegion_, tgtRegion_);

    if (singleMeshProc_ == -1)
    {
        // create global indexing for src and tgt meshes
        globalIndex globalSrcCells(srcRegion_.nCells());
        globalIndex globalTgtCells(tgtRegion_.nCells());

        // Create processor map of overlapping cells. This map gets
        // (possibly remote) cells from the tgt mesh such that they (together)
        // cover all of the src mesh
        autoPtr<mapDistribute> mapPtr = calcProcMap(srcRegion_, tgtRegion_);
        const mapDistribute& map = mapPtr();

        pointField newTgtPoints;
        faceList newTgtFaces;
        labelList newTgtFaceOwners;
        labelList newTgtFaceNeighbours;
        labelList newTgtCellIDs;

        distributeAndMergeCells
        (
            map,
            tgtRegion_,
            globalTgtCells,
            newTgtPoints,
            newTgtFaces,
            newTgtFaceOwners,
            newTgtFaceNeighbours,
            newTgtCellIDs
        );


        // create a new target mesh
        polyMesh newTgt
        (
            IOobject
            (
                "newTgt." + Foam::name(Pstream::myProcNo()),
                tgtRegion_.time().timeName(),
                tgtRegion_.time(),
                IOobject::NO_READ
            ),
            move(newTgtPoints),
            move(newTgtFaces),
            move(newTgtFaceOwners),
            move(newTgtFaceNeighbours),
            false                                   // no parallel comms
        );

        // create some dummy patch info
        List<polyPatch*> patches(1);
        patches[0] = new polyPatch
        (
            "defaultFaces",
            newTgt.nFaces() - newTgt.nInternalFaces(),
            newTgt.nInternalFaces(),
            0,
            newTgt.boundaryMesh(),
            word::null
        );

        newTgt.addPatches(patches);

        // force calculation of tet-base points used for point-in-cell
        (void)newTgt.tetBasePtIs();

        // force construction of cell tree
//        (void)newTgt.cellTree();

        if (debug)
        {
            Pout<< "Created newTgt mesh:" << nl
                << " old cells = " << tgtRegion_.nCells()
                << ", new cells = " << newTgt.nCells() << nl
                << " old faces = " << tgtRegion_.nFaces()
                << ", new faces = " << newTgt.nFaces() << endl;

            if (debug > 1)
            {
                Pout<< "Writing newTgt mesh: " << newTgt.name() << endl;
                newTgt.write();
            }
        }

        calcAddressing(methodName, srcRegion_, newTgt);

        // per source cell the target cell address in newTgt mesh
        forAll(srcToTgtCellAddr_, i)
        {
            labelList& addressing = srcToTgtCellAddr_[i];
            forAll(addressing, addrI)
            {
                addressing[addrI] = newTgtCellIDs[addressing[addrI]];
            }
        }

        // convert target addresses in newTgtMesh into global cell numbering
        forAll(tgtToSrcCellAddr_, i)
        {
            labelList& addressing = tgtToSrcCellAddr_[i];
            forAll(addressing, addrI)
            {
                addressing[addrI] = globalSrcCells.toGlobal(addressing[addrI]);
            }
        }

        // set up as a reverse distribute
        mapDistributeBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            tgtRegion_.nCells(),
            map.constructMap(),
            false,
            map.subMap(),
            false,
            tgtToSrcCellAddr_,
            ListPlusEqOp<label>(),
            flipOp(),
            labelList()
        );

        // set up as a reverse distribute
        mapDistributeBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            tgtRegion_.nCells(),
            map.constructMap(),
            false,
            map.subMap(),
            false,
            tgtToSrcCellWght_,
            ListPlusEqOp<scalar>(),
            flipOp(),
            scalarList()
        );

        // weights normalisation
        normaliseWeights
        (
            "source",
            srcToTgtCellAddr_,
            srcToTgtCellWght_
        );

        normaliseWeights
        (
            "target",
            tgtToSrcCellAddr_,
            tgtToSrcCellWght_
        );

        // cache maps and reset addresses
        List<Map<label>> cMap;
        srcMapPtr_.reset
        (
            new mapDistribute(globalSrcCells, tgtToSrcCellAddr_, cMap)
        );
        tgtMapPtr_.reset
        (
            new mapDistribute(globalTgtCells, srcToTgtCellAddr_, cMap)
        );

        // collect volume intersection contributions
        reduce(V_, sumOp<scalar>());
    }
    else
    {
        calcAddressing(methodName, srcRegion_, tgtRegion_);

        normaliseWeights
        (
            "source",
            srcToTgtCellAddr_,
            srcToTgtCellWght_
        );

        normaliseWeights
        (
            "target",
            tgtToSrcCellAddr_,
            tgtToSrcCellWght_
        );
    }

    Info<< "    Overlap volume: " << V_ << endl;
}


Foam::AMIInterpolation::interpolationMethod
Foam::meshToMesh::interpolationMethodAMI(const interpolationMethod method)
{
    switch (method)
    {
        case imDirect:
        {
            return AMIInterpolation::imDirect;
            break;
        }
        case imMapNearest:
        {
            return AMIInterpolation::imMapNearest;
            break;
        }
        case imCellVolumeWeight:
        {
            return AMIInterpolation::imFaceAreaWeight;
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration " << method
                << abort(FatalError);
        }
    }

    return AMIInterpolation::imDirect;
}


void Foam::meshToMesh::calculatePatchAMIs(const word& AMIMethodName)
{
    if (!patchAMIs_.empty())
    {
        FatalErrorInFunction
            << "patch AMI already calculated"
            << exit(FatalError);
    }

    patchAMIs_.setSize(srcPatchID_.size());

    forAll(srcPatchID_, i)
    {
        label srcPatchi = srcPatchID_[i];
        label tgtPatchi = tgtPatchID_[i];

        const polyPatch& srcPP = srcRegion_.boundaryMesh()[srcPatchi];
        const polyPatch& tgtPP = tgtRegion_.boundaryMesh()[tgtPatchi];

        Info<< "Creating AMI between source patch " << srcPP.name()
            << " and target patch " << tgtPP.name()
            << " using " << AMIMethodName
            << endl;

        Info<< incrIndent;

        patchAMIs_.set
        (
            i,
            new AMIInterpolation
            (
                srcPP,
                tgtPP,
                faceAreaIntersect::tmMesh,
                false,
                AMIMethodName,
                -1,
                true // flip target patch since patch normals are aligned
            )
        );

        Info<< decrIndent;
    }
}


void Foam::meshToMesh::constructNoCuttingPatches
(
    const word& methodName,
    const word& AMIMethodName,
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
    calculatePatchAMIs(AMIMethodName);
}


void Foam::meshToMesh::constructFromCuttingPatches
(
    const word& methodName,
    const word& AMIMethodName,
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
    calculatePatchAMIs(AMIMethodName);

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
    const interpolationMethod& method,
    bool interpAllPatches
)
:
    srcRegion_(src),
    tgtRegion_(tgt),
    srcPatchID_(),
    tgtPatchID_(),
    patchAMIs_(),
    cuttingPatches_(),
    srcToTgtCellAddr_(),
    tgtToSrcCellAddr_(),
    srcToTgtCellWght_(),
    tgtToSrcCellWght_(),
    V_(0.0),
    singleMeshProc_(-1),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr)
{
    constructNoCuttingPatches
    (
        interpolationMethodNames_[method],
        AMIInterpolation::interpolationMethodToWord
        (
            interpolationMethodAMI(method)
        ),
        interpAllPatches
    );
}


Foam::meshToMesh::meshToMesh
(
    const polyMesh& src,
    const polyMesh& tgt,
    const word& methodName,
    const word& AMIMethodName,
    bool interpAllPatches
)
:
    srcRegion_(src),
    tgtRegion_(tgt),
    srcPatchID_(),
    tgtPatchID_(),
    patchAMIs_(),
    cuttingPatches_(),
    srcToTgtCellAddr_(),
    tgtToSrcCellAddr_(),
    srcToTgtCellWght_(),
    tgtToSrcCellWght_(),
    V_(0.0),
    singleMeshProc_(-1),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr)
{
    constructNoCuttingPatches(methodName, AMIMethodName, interpAllPatches);
}


Foam::meshToMesh::meshToMesh
(
    const polyMesh& src,
    const polyMesh& tgt,
    const interpolationMethod& method,
    const HashTable<word>& patchMap,
    const wordList& cuttingPatches
)
:
    srcRegion_(src),
    tgtRegion_(tgt),
    srcPatchID_(),
    tgtPatchID_(),
    patchAMIs_(),
    cuttingPatches_(),
    srcToTgtCellAddr_(),
    tgtToSrcCellAddr_(),
    srcToTgtCellWght_(),
    tgtToSrcCellWght_(),
    V_(0.0),
    singleMeshProc_(-1),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr)
{
    constructFromCuttingPatches
    (
        interpolationMethodNames_[method],
        AMIInterpolation::interpolationMethodToWord
        (
            interpolationMethodAMI(method)
        ),
        patchMap,
        cuttingPatches
    );
}


Foam::meshToMesh::meshToMesh
(
    const polyMesh& src,
    const polyMesh& tgt,
    const word& methodName,     // internal mapping
    const word& AMIMethodName,  // boundary mapping
    const HashTable<word>& patchMap,
    const wordList& cuttingPatches
)
:
    srcRegion_(src),
    tgtRegion_(tgt),
    srcPatchID_(),
    tgtPatchID_(),
    patchAMIs_(),
    cuttingPatches_(),
    srcToTgtCellAddr_(),
    tgtToSrcCellAddr_(),
    srcToTgtCellWght_(),
    tgtToSrcCellWght_(),
    V_(0.0),
    singleMeshProc_(-1),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr)
{
    constructFromCuttingPatches
    (
        methodName,
        AMIMethodName,
        patchMap,
        cuttingPatches
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshToMesh::~meshToMesh()
{}


// ************************************************************************* //
