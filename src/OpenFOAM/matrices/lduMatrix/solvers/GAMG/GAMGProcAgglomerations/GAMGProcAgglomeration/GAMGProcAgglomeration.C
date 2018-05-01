/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "GAMGProcAgglomeration.H"
#include "GAMGAgglomeration.H"
#include "lduMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GAMGProcAgglomeration, 0);
    defineRunTimeSelectionTable(GAMGProcAgglomeration, GAMGAgglomeration);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::GAMGProcAgglomeration::printStats
(
    Ostream& os,
    GAMGAgglomeration& agglom
) const
{
    for (label levelI = 0; levelI <= agglom.size(); levelI++)
    {
        if (agglom.hasMeshLevel(levelI))
        {
            os  << agglom.meshLevel(levelI).info() << endl;
        }
        else
        {
            os  << "Level " << levelI << " has no fine mesh:" << endl;
        }

        if
        (
            levelI < agglom.restrictAddressing_.size()
         && agglom.restrictAddressing_.set(levelI)
        )
        {
            const labelList& cellRestrict =
                agglom.restrictAddressing(levelI);
            const labelList& faceRestrict =
                agglom.faceRestrictAddressing(levelI);

            os  << "Level " << levelI << " agglomeration:" << nl
                << "    nCoarseCells:" << agglom.nCells(levelI) << nl
                << "    nCoarseFaces:" << agglom.nFaces(levelI) << nl
                << "    cellRestriction:"
                << " size:" << cellRestrict.size()
                << " max:" << max(cellRestrict)
                << nl
                << "    faceRestriction:"
                << " size:" << faceRestrict.size()
                << " max:" << max(faceRestrict)
                << nl;

            const labelListList& patchFaceRestrict =
                agglom.patchFaceRestrictAddressing(levelI);
            forAll(patchFaceRestrict, i)
            {
                if (patchFaceRestrict[i].size())
                {
                    const labelList& faceRestrict =
                        patchFaceRestrict[i];
                    os  << "        " << i
                        << " size:" << faceRestrict.size()
                        << " max:" << max(faceRestrict)
                        << nl;
                }
            }
        }
        if
        (
            levelI < agglom.procCellOffsets_.size()
         && agglom.procCellOffsets_.set(levelI)
        )
        {
            os  << "    procCellOffsets:" << agglom.procCellOffsets_[levelI]
                << nl
                << "    procAgglomMap:" << agglom.procAgglomMap_[levelI]
                << nl
                << "    procIDs:" << agglom.agglomProcIDs_[levelI]
                << nl
                << "    comm:" << agglom.procCommunicator_[levelI]
                << endl;
        }

        os  << endl;
    }
    os  << endl;
}


Foam::labelListList Foam::GAMGProcAgglomeration::globalCellCells
(
    const lduMesh& mesh
)
{
    const lduAddressing& addr = mesh.lduAddr();
    lduInterfacePtrsList interfaces = mesh.interfaces();

    const label myProcID = Pstream::myProcNo(mesh.comm());

    globalIndex globalNumbering
    (
        addr.size(),
        Pstream::msgType(),
        mesh.comm(),
        Pstream::parRun()
    );

    labelList globalIndices(addr.size());
    forAll(globalIndices, celli)
    {
        globalIndices[celli] = globalNumbering.toGlobal(myProcID, celli);
    }


    // Get the interface cells
    PtrList<labelList> nbrGlobalCells(interfaces.size());
    {
        // Initialise transfer of restrict addressing on the interface
        forAll(interfaces, inti)
        {
            if (interfaces.set(inti))
            {
                interfaces[inti].initInternalFieldTransfer
                (
                    Pstream::commsTypes::nonBlocking,
                    globalIndices
                );
            }
        }

        if (Pstream::parRun())
        {
            Pstream::waitRequests();
        }

        forAll(interfaces, inti)
        {
            if (interfaces.set(inti))
            {
                nbrGlobalCells.set
                (
                    inti,
                    new labelList
                    (
                        interfaces[inti].internalFieldTransfer
                        (
                            Pstream::commsTypes::nonBlocking,
                            globalIndices
                        )
                    )
                );
            }
        }
    }


    // Scan the neighbour list to find out how many times the cell
    // appears as a neighbour of the face. Done this way to avoid guessing
    // and resizing list
    labelList nNbrs(addr.size(), 1);

    const labelUList& nbr = addr.upperAddr();
    const labelUList& own = addr.lowerAddr();

    {
        forAll(nbr, facei)
        {
            nNbrs[nbr[facei]]++;
            nNbrs[own[facei]]++;
        }

        forAll(interfaces, inti)
        {
            if (interfaces.set(inti))
            {
                const labelUList& faceCells = interfaces[inti].faceCells();

                forAll(faceCells, i)
                {
                    nNbrs[faceCells[i]]++;
                }
            }
        }
    }


    // Create cell-cells addressing
    labelListList cellCells(addr.size());

    forAll(cellCells, celli)
    {
        cellCells[celli].setSize(nNbrs[celli], -1);
    }

    // Reset the list of number of neighbours to zero
    nNbrs = 0;

    // Scatter the neighbour faces
    forAll(nbr, facei)
    {
        label c0 = own[facei];
        label c1 = nbr[facei];

        cellCells[c0][nNbrs[c0]++] = globalIndices[c1];
        cellCells[c1][nNbrs[c1]++] = globalIndices[c0];
    }
    forAll(interfaces, inti)
    {
        if (interfaces.set(inti))
        {
            const labelUList& faceCells = interfaces[inti].faceCells();

            forAll(faceCells, i)
            {
                label c0 = faceCells[i];
                cellCells[c0][nNbrs[c0]++] = nbrGlobalCells[inti][i];
            }
        }
    }

    forAll(cellCells, celli)
    {
        Foam::stableSort(cellCells[celli]);
    }

    // Replace the initial element (always -1) with the local cell
    forAll(cellCells, celli)
    {
        cellCells[celli][0] = globalIndices[celli];
    }

    return cellCells;
}


bool Foam::GAMGProcAgglomeration::agglomerate
(
    const label fineLevelIndex,
    const labelList& procAgglomMap,
    const labelList& masterProcs,
    const List<label>& agglomProcIDs,
    const label procAgglomComm
)
{
    const lduMesh& levelMesh = agglom_.meshLevels_[fineLevelIndex];
    label levelComm = levelMesh.comm();

    if (Pstream::myProcNo(levelComm) != -1)
    {
        // Collect meshes and restrictAddressing onto master
        // Overwrites the fine mesh (meshLevels_[index-1]) and addressing
        // from fine mesh to coarse mesh (restrictAddressing_[index]).
        agglom_.procAgglomerateLduAddressing
        (
            levelComm,
            procAgglomMap,
            agglomProcIDs,
            procAgglomComm,

            fineLevelIndex               // fine level index
        );

        // Combine restrict addressing only onto master
        for
        (
            label levelI = fineLevelIndex+1;
            levelI < agglom_.meshLevels_.size();
            levelI++
        )
        {
            agglom_.procAgglomerateRestrictAddressing
            (
                levelComm,
                agglomProcIDs,
                levelI
            );
        }

        if (Pstream::myProcNo(levelComm) == agglomProcIDs[0])
        {
            // On master. Recreate coarse meshes from restrict addressing
            for
            (
                label levelI = fineLevelIndex;
                levelI < agglom_.meshLevels_.size();
                levelI++
            )
            {
                agglom_.agglomerateLduAddressing(levelI);
            }
        }
        else
        {
            // Agglomerated away. Clear mesh storage.
            for
            (
                label levelI = fineLevelIndex+1;
                levelI <= agglom_.size();
                levelI++
            )
            {
                agglom_.clearLevel(levelI);
            }
        }
    }

    // Should check!
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GAMGProcAgglomeration::GAMGProcAgglomeration
(
    GAMGAgglomeration& agglom,
    const dictionary& controlDict
)
:
    agglom_(agglom)
{}


Foam::autoPtr<Foam::GAMGProcAgglomeration> Foam::GAMGProcAgglomeration::New
(
    const word& type,
    GAMGAgglomeration& agglom,
    const dictionary& controlDict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing GAMGProcAgglomeration" << endl;
    }

    GAMGAgglomerationConstructorTable::iterator cstrIter =
        GAMGAgglomerationConstructorTablePtr_->find(type);

    if (cstrIter == GAMGAgglomerationConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown GAMGProcAgglomeration type "
            << type << " for GAMGAgglomeration " << agglom.type() << nl << nl
            << "Valid GAMGProcAgglomeration types are :" << endl
            << GAMGAgglomerationConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<GAMGProcAgglomeration>(cstrIter()(agglom, controlDict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GAMGProcAgglomeration::~GAMGProcAgglomeration()
{}


// ************************************************************************* //
