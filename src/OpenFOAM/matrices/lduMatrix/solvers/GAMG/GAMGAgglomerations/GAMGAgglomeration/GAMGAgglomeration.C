/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "GAMGAgglomeration.H"
#include "lduMesh.H"
#include "lduMatrix.H"
#include "Time.H"
#include "GAMGInterface.H"
#include "GAMGProcAgglomeration.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GAMGAgglomeration, 0);
    defineRunTimeSelectionTable(GAMGAgglomeration, lduMesh);
    defineRunTimeSelectionTable(GAMGAgglomeration, lduMatrix);
    defineRunTimeSelectionTable(GAMGAgglomeration, geometry);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GAMGAgglomeration::compactLevels(const label nCreatedLevels)
{
    nCells_.setSize(nCreatedLevels);
    restrictAddressing_.setSize(nCreatedLevels);
    nFaces_.setSize(nCreatedLevels);
    faceRestrictAddressing_.setSize(nCreatedLevels);
    faceFlipMap_.setSize(nCreatedLevels);
    nPatchFaces_.setSize(nCreatedLevels);
    patchFaceRestrictAddressing_.setSize(nCreatedLevels);
    meshLevels_.setSize(nCreatedLevels);

    // Have procCommunicator_ always, even if not procAgglomerating
    procCommunicator_.setSize(nCreatedLevels + 1);
    if (processorAgglomerate())
    {
        procAgglomMap_.setSize(nCreatedLevels);
        agglomProcIDs_.setSize(nCreatedLevels);
        procCellOffsets_.setSize(nCreatedLevels);
        procFaceMap_.setSize(nCreatedLevels);
        procBoundaryMap_.setSize(nCreatedLevels);
        procBoundaryFaceMap_.setSize(nCreatedLevels);

        procAgglomeratorPtr_().agglomerate();


    }

    // Print a bit
    if (processorAgglomerate() && debug)
    {
        Info<< "GAMGAgglomeration:" << nl
            << "    local agglomerator     : " << type() << nl;
        if (processorAgglomerate())
        {
            Info<< "    processor agglomerator : "
                << procAgglomeratorPtr_().type() << nl
                << nl;
        }

        Info<< setw(36) << "nCells"
            << setw(20) << "nFaces/nCells"
            << setw(20) << "nInterfaces"
            << setw(20) << "nIntFaces/nCells"
            << setw(12) << "profile"
            << nl
            << setw(8) << "Level"
            << setw(8) << "nProcs"
            << "    "
            << setw(8) << "avg"
            << setw(8) << "max"
            << "    "
            << setw(8) << "avg"
            << setw(8) << "max"
            << "    "
            << setw(8) << "avg"
            << setw(8) << "max"
            << "    "
            << setw(8) << "avg"
            << setw(8) << "max"
            //<< "    "
            << setw(12) << "avg"
            << nl
            << setw(8) << "-----"
            << setw(8) << "------"
            << "    "
            << setw(8) << "---"
            << setw(8) << "---"
            << "    "
            << setw(8) << "---"
            << setw(8) << "---"
            << "    "
            << setw(8) << "---"
            << setw(8) << "---"
            << "    "
            << setw(8) << "---"
            << setw(8) << "---"
            //<< "    "
            << setw(12) << "---"
            //<< "    "
            << nl;

        for (label levelI = 0; levelI <= size(); levelI++)
        {
            label nProcs = 0;
            label nCells = 0;
            scalar faceCellRatio = 0;
            label nInterfaces = 0;
            label nIntFaces = 0;
            scalar ratio = 0.0;
            scalar profile = 0.0;

            if (hasMeshLevel(levelI))
            {
                nProcs = 1;

                const lduMesh& fineMesh = meshLevel(levelI);
                nCells = fineMesh.lduAddr().size();
                faceCellRatio =
                    scalar(fineMesh.lduAddr().lowerAddr().size())/nCells;

                const lduInterfacePtrsList interfaces =
                    fineMesh.interfaces();
                forAll(interfaces, i)
                {
                    if (interfaces.set(i))
                    {
                        nInterfaces++;
                        nIntFaces += interfaces[i].faceCells().size();
                    }
                }
                ratio = scalar(nIntFaces)/nCells;

                profile = fineMesh.lduAddr().band().second();
            }

            label totNprocs = returnReduce(nProcs, sumOp<label>());

            label maxNCells = returnReduce(nCells, maxOp<label>());
            label totNCells = returnReduce(nCells, sumOp<label>());

            scalar maxFaceCellRatio =
                returnReduce(faceCellRatio, maxOp<scalar>());
            scalar totFaceCellRatio =
                returnReduce(faceCellRatio, sumOp<scalar>());

            label maxNInt = returnReduce(nInterfaces, maxOp<label>());
            label totNInt = returnReduce(nInterfaces, sumOp<label>());

            scalar maxRatio = returnReduce(ratio, maxOp<scalar>());
            scalar totRatio = returnReduce(ratio, sumOp<scalar>());

            scalar totProfile = returnReduce(profile, sumOp<scalar>());

            int oldPrecision = Info().precision(4);

            Info<< setw(8) << levelI
                << setw(8) << totNprocs
                << "    "
                << setw(8) << totNCells/totNprocs
                << setw(8) << maxNCells
                << "    "
                << setw(8) << totFaceCellRatio/totNprocs
                << setw(8) << maxFaceCellRatio
                << "    "
                << setw(8) << scalar(totNInt)/totNprocs
                << setw(8) << maxNInt
                << "    "
                << setw(8) << totRatio/totNprocs
                << setw(8) << maxRatio
                << setw(12) << totProfile/totNprocs
                << nl;

            Info().precision(oldPrecision);
        }
        Info<< endl;
    }
}


bool Foam::GAMGAgglomeration::continueAgglomerating
(
    const label nCoarseCells
) const
{
    // Check the need for further agglomeration on all processors
    bool contAgg = nCoarseCells >= nCellsInCoarsestLevel_;
    mesh().reduce(contAgg, andOp<bool>());
    return contAgg;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GAMGAgglomeration::GAMGAgglomeration
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
:
    MeshObject<lduMesh, Foam::GeometricMeshObject, GAMGAgglomeration>(mesh),

    maxLevels_(50),

    nCellsInCoarsestLevel_
    (
        readLabel(controlDict.lookup("nCellsInCoarsestLevel"))
    ),
    meshInterfaces_(mesh.interfaces()),
    procAgglomeratorPtr_
    (
        (
            (UPstream::nProcs(mesh.comm()) > 1)
         && controlDict.found("processorAgglomerator")
        )
      ? GAMGProcAgglomeration::New
        (
            controlDict.lookup("processorAgglomerator"),
            *this,
            controlDict
        )
      : autoPtr<GAMGProcAgglomeration>(NULL)
    ),

    nCells_(maxLevels_),
    restrictAddressing_(maxLevels_),
    nFaces_(maxLevels_),
    faceRestrictAddressing_(maxLevels_),
    faceFlipMap_(maxLevels_),
    nPatchFaces_(maxLevels_),
    patchFaceRestrictAddressing_(maxLevels_),

    meshLevels_(maxLevels_)
{
    procCommunicator_.setSize(maxLevels_ + 1, -1);
    if (processorAgglomerate())
    {
        procAgglomMap_.setSize(maxLevels_);
        agglomProcIDs_.setSize(maxLevels_);
        procCellOffsets_.setSize(maxLevels_);
        procFaceMap_.setSize(maxLevels_);
        procBoundaryMap_.setSize(maxLevels_);
        procBoundaryFaceMap_.setSize(maxLevels_);
    }
}


const Foam::GAMGAgglomeration& Foam::GAMGAgglomeration::New
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
{
    if
    (
        !mesh.thisDb().foundObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        )
    )
    {
        const word agglomeratorType(controlDict.lookup("agglomerator"));

        const_cast<Time&>(mesh.thisDb().time()).libs().open
        (
            controlDict,
            "geometricGAMGAgglomerationLibs",
            lduMeshConstructorTablePtr_
        );

        lduMeshConstructorTable::iterator cstrIter =
            lduMeshConstructorTablePtr_->find(agglomeratorType);

        if (cstrIter == lduMeshConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "GAMGAgglomeration::New"
                "(const lduMesh& mesh, const dictionary& controlDict)"
            )   << "Unknown GAMGAgglomeration type "
                << agglomeratorType << ".\n"
                << "Valid matrix GAMGAgglomeration types are :"
                << lduMatrixConstructorTablePtr_->sortedToc() << endl
                << "Valid geometric GAMGAgglomeration types are :"
                << lduMeshConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return store(cstrIter()(mesh, controlDict).ptr());
    }
    else
    {
        return mesh.thisDb().lookupObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        );
    }
}


const Foam::GAMGAgglomeration& Foam::GAMGAgglomeration::New
(
    const lduMatrix& matrix,
    const dictionary& controlDict
)
{
    const lduMesh& mesh = matrix.mesh();

    if
    (
        !mesh.thisDb().foundObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        )
    )
    {
        const word agglomeratorType(controlDict.lookup("agglomerator"));

        const_cast<Time&>(mesh.thisDb().time()).libs().open
        (
            controlDict,
            "algebraicGAMGAgglomerationLibs",
            lduMatrixConstructorTablePtr_
        );

        if
        (
            !lduMatrixConstructorTablePtr_
         || !lduMatrixConstructorTablePtr_->found(agglomeratorType)
        )
        {
            return New(mesh, controlDict);
        }
        else
        {
            lduMatrixConstructorTable::iterator cstrIter =
                lduMatrixConstructorTablePtr_->find(agglomeratorType);

            return store(cstrIter()(matrix, controlDict).ptr());
        }
    }
    else
    {
        return mesh.thisDb().lookupObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        );
    }
}


Foam::autoPtr<Foam::GAMGAgglomeration> Foam::GAMGAgglomeration::New
(
    const lduMesh& mesh,
    const scalarField& cellVolumes,
    const vectorField& faceAreas,
    const dictionary& controlDict
)
{
    const word agglomeratorType(controlDict.lookup("agglomerator"));

    const_cast<Time&>(mesh.thisDb().time()).libs().open
    (
        controlDict,
        "geometricGAMGAgglomerationLibs",
        geometryConstructorTablePtr_
    );

    geometryConstructorTable::iterator cstrIter =
        geometryConstructorTablePtr_->find(agglomeratorType);

    if (cstrIter == geometryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "GAMGAgglomeration::New"
            "(const lduMesh& mesh, const scalarField&"
            ", const vectorField&, const dictionary& controlDict)"
        )   << "Unknown GAMGAgglomeration type "
            << agglomeratorType << ".\n"
            << "Valid geometric GAMGAgglomeration types are :"
            << geometryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<GAMGAgglomeration>
    (
        cstrIter()
        (
            mesh,
            cellVolumes,
            faceAreas,
            controlDict
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GAMGAgglomeration::~GAMGAgglomeration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::lduMesh& Foam::GAMGAgglomeration::meshLevel
(
    const label i
) const
{
    if (i == 0)
    {
        return mesh_;
    }
    else
    {
        return meshLevels_[i - 1];
    }
}


bool Foam::GAMGAgglomeration::hasMeshLevel(const label i) const
{
    if (i == 0)
    {
        return true;
    }
    else
    {
        return meshLevels_.set(i - 1);
    }
}


const Foam::lduInterfacePtrsList& Foam::GAMGAgglomeration::interfaceLevel
(
    const label i
) const
{
    if (i == 0)
    {
        return meshInterfaces_;
    }
    else
    {
        return meshLevels_[i - 1].rawInterfaces();
    }
}


void Foam::GAMGAgglomeration::clearLevel(const label i)
{
    if (hasMeshLevel(i))
    {
        meshLevels_.set(i - 1, NULL);

        if (i < nCells_.size())
        {
            nCells_[i] = -555;
            restrictAddressing_.set(i, NULL);
            nFaces_[i] = -666;
            faceRestrictAddressing_.set(i, NULL);
            faceFlipMap_.set(i, NULL);
            nPatchFaces_.set(i, NULL);
            patchFaceRestrictAddressing_.set(i, NULL);
        }
    }
}


const Foam::labelList& Foam::GAMGAgglomeration::procAgglomMap
(
    const label leveli
) const
{
    return procAgglomMap_[leveli];
}


const Foam::labelList& Foam::GAMGAgglomeration::agglomProcIDs
(
    const label leveli
) const
{
    return agglomProcIDs_[leveli];
}


bool Foam::GAMGAgglomeration::hasProcMesh(const label leveli) const
{
    return procCommunicator_[leveli] != -1;
}


Foam::label Foam::GAMGAgglomeration::procCommunicator(const label leveli) const
{
    return procCommunicator_[leveli];
}


const Foam::labelList& Foam::GAMGAgglomeration::cellOffsets
(
    const label leveli
) const
{
    return procCellOffsets_[leveli];
}


const Foam::labelListList& Foam::GAMGAgglomeration::faceMap
(
    const label leveli
) const
{
    return procFaceMap_[leveli];
}


const Foam::labelListList& Foam::GAMGAgglomeration::boundaryMap
(
    const label leveli
) const
{
    return procBoundaryMap_[leveli];
}


const Foam::labelListListList& Foam::GAMGAgglomeration::boundaryFaceMap
(
    const label leveli
) const
{
    return procBoundaryFaceMap_[leveli];
}


bool Foam::GAMGAgglomeration::checkRestriction
(
    labelList& newRestrict,
    label& nNewCoarse,
    const lduAddressing& fineAddressing,
    const labelUList& restrict,
    const label nCoarse
)
{
    if (fineAddressing.size() != restrict.size())
    {
        FatalErrorIn
        (
            "checkRestriction(..)"
        )   << "nCells:" << fineAddressing.size()
            << " agglom:" << restrict.size()
            << abort(FatalError);
    }

    // Seed (master) for every region
    labelList master(identity(fineAddressing.size()));

    // Now loop and transport master through region
    const labelUList& lower = fineAddressing.lowerAddr();
    const labelUList& upper = fineAddressing.upperAddr();

    while (true)
    {
        label nChanged = 0;

        forAll(lower, faceI)
        {
            label own = lower[faceI];
            label nei = upper[faceI];

            if (restrict[own] == restrict[nei])
            {
                // coarse-mesh-internal face

                if (master[own] < master[nei])
                {
                    master[nei] = master[own];
                    nChanged++;
                }
                else if (master[own] > master[nei])
                {
                    master[own] = master[nei];
                    nChanged++;
                }
            }
        }

        reduce(nChanged, sumOp<label>());

        if (nChanged == 0)
        {
            break;
        }
    }


    // Count number of regions/masters per coarse cell
    labelListList coarseToMasters(nCoarse);
    nNewCoarse = 0;
    forAll(restrict, cellI)
    {
        labelList& masters = coarseToMasters[restrict[cellI]];

        if (findIndex(masters, master[cellI]) == -1)
        {
            masters.append(master[cellI]);
            nNewCoarse++;
        }
    }

    if (nNewCoarse > nCoarse)
    {
        //WarningIn("GAMGAgglomeration::checkRestriction(..)")
        //    << "Have " << nCoarse
        //    << " agglomerated cells but " << nNewCoarse
        //    << " disconnected regions" << endl;

        // Keep coarseToMasters[0] the original coarse, allocate new ones
        // for the others
        labelListList coarseToNewCoarse(coarseToMasters.size());

        nNewCoarse = nCoarse;

        forAll(coarseToMasters, coarseI)
        {
            const labelList& masters = coarseToMasters[coarseI];

            labelList& newCoarse = coarseToNewCoarse[coarseI];
            newCoarse.setSize(masters.size());
            newCoarse[0] = coarseI;
            for (label i = 1; i < newCoarse.size(); i++)
            {
                newCoarse[i] = nNewCoarse++;
            }
        }

        newRestrict.setSize(fineAddressing.size());
        forAll(restrict, cellI)
        {
            label coarseI = restrict[cellI];

            label index = findIndex(coarseToMasters[coarseI], master[cellI]);
            newRestrict[cellI] = coarseToNewCoarse[coarseI][index];
        }

        return false;
    }
    else
    {
        return true;
    }
}


// ************************************************************************* //
