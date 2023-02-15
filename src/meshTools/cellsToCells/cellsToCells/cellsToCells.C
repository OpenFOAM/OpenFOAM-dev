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

#include "cellsToCells.H"
#include "globalIndex.H"
#include "PatchTools.H"
#include "patchToPatchTools.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "processorPolyPatch.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellsToCells, 0);
    defineRunTimeSelectionTable(cellsToCells, word);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::cellsToCells::initialise
(
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh
)
{
    srcLocalTgtCells_.setSize(srcMesh.nCells());
    srcWeights_.setSize(srcMesh.nCells());
    forAll(srcLocalTgtCells_, srcCelli)
    {
        srcLocalTgtCells_[srcCelli].clear();
        srcWeights_[srcCelli].clear();
    }

    tgtLocalSrcCells_.setSize(tgtMesh.nCells());
    tgtWeights_.setSize(tgtMesh.nCells());
    forAll(tgtLocalSrcCells_, tgtCelli)
    {
        tgtLocalSrcCells_[tgtCelli].clear();
        tgtWeights_[tgtCelli].clear();
    }
}


Foam::labelList Foam::cellsToCells::maskCells
(
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh
) const
{
    boundBox meshBb
    (
        max(srcMesh.bounds().min(), tgtMesh.bounds().min()),
        min(srcMesh.bounds().max(), tgtMesh.bounds().max())
    );

    meshBb.inflate(0.01);

    const cellList& srcCells = srcMesh.cells();
    const faceList& srcFaces = srcMesh.faces();
    const pointField& srcPts = srcMesh.points();

    DynamicList<label> resultDyn(srcMesh.nCells());
    forAll(srcCells, srcCelli)
    {
        const boundBox cellBb
        (
            srcCells[srcCelli].points(srcFaces, srcPts),
            false
        );

        if (meshBb.overlaps(cellBb))
        {
            resultDyn.append(srcCelli);
        }
    }

    labelList result;
    result.transfer(resultDyn);
    return result;
}


void Foam::cellsToCells::appendNbrCells
(
    const label celli,
    const polyMesh& mesh,
    const DynamicList<label>& visitedCells,
    DynamicList<label>& nbrCells
) const
{
    // Get all cell-cells
    const labelList& allNbrCells = mesh.cellCells()[celli];

    // Filter out cells already visited
    forAll(allNbrCells, i)
    {
        const label nbrCelli = allNbrCells[i];

        if
        (
            findIndex(visitedCells, nbrCelli) == -1
         && findIndex(nbrCells, nbrCelli) == -1
        )
        {
            nbrCells.append(nbrCelli);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellsToCells::cellsToCells()
:
    singleProcess_(-1),
    srcLocalTgtCells_(),
    tgtLocalSrcCells_(),
    srcWeights_(),
    tgtWeights_(),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr),
    localSrcProcCellsPtr_(nullptr),
    localTgtProcCellsPtr_(nullptr),
    localTgtMeshPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellsToCells::~cellsToCells()
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::cellsToCells> Foam::cellsToCells::New
(
    const word& cellsToCellsType
)
{
    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(cellsToCellsType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName << " type "
            << cellsToCellsType << endl << endl
            << "Valid " << typeName << " types are : " << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::PackedBoolList Foam::cellsToCells::srcCoupled() const
{
    PackedBoolList result(srcLocalTgtCells_.size());
    forAll(srcLocalTgtCells_, srcCelli)
    {
        result[srcCelli] = !srcLocalTgtCells_[srcCelli].empty();
    }
    return result;
}


Foam::PackedBoolList Foam::cellsToCells::tgtCoupled() const
{
    PackedBoolList result(tgtLocalSrcCells_.size());
    forAll(tgtLocalSrcCells_, tgtCelli)
    {
        result[tgtCelli] = !tgtLocalSrcCells_[tgtCelli].empty();
    }
    return result;
}


Foam::remote Foam::cellsToCells::srcToTgtPoint
(
    const polyMesh& tgtMesh,
    const label srcCelli,
    const point& p
) const
{
    forAll(srcLocalTgtCells_[srcCelli], i)
    {
        const label tgtCelli = srcLocalTgtCells_[srcCelli][i];

        const polyMesh& localTgtMesh =
            singleProcess_ == -1 ? localTgtMeshPtr_() : tgtMesh;

        if (localTgtMesh.pointInCell(p, tgtCelli))
        {
            return
                singleProcess_ == -1
              ? localTgtProcCellsPtr_()[tgtCelli]
              : remote(Pstream::myProcNo(), tgtCelli);
        }
    }

    return remote();
}


Foam::scalar Foam::cellsToCells::update
(
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh
)
{
    cpuTime time;

    // Determine numbers of faces on both sides, report, and quit if either
    // side is empty
    const label srcTotalSize = returnReduce(srcMesh.nCells(), sumOp<label>());
    const label tgtTotalSize = returnReduce(tgtMesh.nCells(), sumOp<label>());
    if (srcTotalSize == 0 || tgtTotalSize == 0)
    {
        return 0;
    }

    Info<< indent << typeName << ": Calculating couplings between "
        << srcTotalSize << " source cells and " << tgtTotalSize
        << " target cells" << incrIndent << endl;

    singleProcess_ =
        patchToPatchTools::singleProcess
        (
            srcMesh.nCells(),
            tgtMesh.nCells()
        );

    scalar V = 0;

    if (isSingleProcess())
    {
        // Do the intersection
        V = calculate(srcMesh, tgtMesh);

        // Normalise the weights
        normalise(srcMesh, srcLocalTgtCells_, srcWeights_);
        normalise(tgtMesh, tgtLocalSrcCells_, tgtWeights_);
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
                tgtMeshSendCells(srcMesh, tgtMesh)
            );
        localTgtProcCellsPtr_.reset
        (
            new List<remote>
            (
                distributeMesh
                (
                    tgtMapPtr_(),
                    tgtMesh,
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
        V = calculate(srcMesh, localTgtMesh);

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
            tgtMesh.nCells(),
            tgtMapPtr_(),
            localSrcProcCellsPtr_(),
            tgtLocalSrcCells_
        );

        // Collect the weights on the target
        patchToPatchTools::rDistributeListList
        (
            tgtMesh.nCells(),
            tgtMapPtr_(),
            tgtWeights_
        );

        // Normalise the weights
        normalise(srcMesh, srcLocalTgtCells_, srcWeights_);
        normalise(tgtMesh, tgtLocalSrcCells_, tgtWeights_);

        // Collect volume intersection contributions
        reduce(V, sumOp<scalar>());
    }

    label nCouples = 0;
    forAll(srcLocalTgtCells_, srcCelli)
    {
        nCouples += srcLocalTgtCells_[srcCelli].size();
    }
    forAll(tgtLocalSrcCells_, tgtCelli)
    {
        nCouples += tgtLocalSrcCells_[tgtCelli].size();
    }
    reduce(nCouples, sumOp<label>());

    if (nCouples != 0)
    {
        Info<< indent << "Overlapping volume = " << V << endl
            << indent << nCouples << " couplings calculated in "
            << time.cpuTimeIncrement() << 's' << endl;
    }
    else
    {
        Info<< indent << "No couplings found" << endl;
    }

    Info<< decrIndent;

    return V;
}


// ************************************************************************* //
