/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2026 OpenFOAM Foundation
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

#include "cellsToCellsExtrapolation.H"
#include "globalIndex.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellsToCellsExtrapolation, 0);
    defineRunTimeSelectionTable(cellsToCellsExtrapolation, word);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::cellsToCellsExtrapolation::update
(
    const PackedBoolList& cellCoupleds
)
{
    uncoupledCells_ =
        selectIndices(cellCoupleds, [](const bool b) { return !b ;});

    extrapolation_ = uncoupledCells_.size();

    reduce(extrapolation_, orOp());
}


void Foam::cellsToCellsExtrapolation::getCutFaces
(
    const polyMesh& mesh,
    const PackedBoolList& cellCoupleds,
    const globalIndex& globalCellIndex,
    labelList& cutFaces,
    labelList& cutFaceCoupledCells
) const
{
    // Get some information regarding the cells on the other side of couplings
    labelList bFaceNbrCells(mesh.nFaces() - mesh.nInternalFaces());
    boolList bFaceNbrIsCoupled(mesh.nFaces() - mesh.nInternalFaces());
    forAll(bFaceNbrIsCoupled, bFacei)
    {
        const label owni = mesh.faceOwner()[bFacei + mesh.nInternalFaces()];
        bFaceNbrCells[bFacei] = globalCellIndex.toGlobal(owni);
        bFaceNbrIsCoupled[bFacei] = cellCoupleds[owni];
    }
    syncTools::swapBoundaryFaceList(mesh, bFaceNbrCells);
    syncTools::swapBoundaryFaceList(mesh, bFaceNbrIsCoupled);

    // Determine the "cut" faces that separate coupled and un-coupled cells
    DynamicList<label> cutFacesDyn, cutFacesCoupledCellDyn;
    for (label facei = 0; facei < mesh.nFaces(); ++ facei)
    {
        const label owni = mesh.faceOwner()[facei];
        const bool ownIsCoupled = cellCoupleds[owni];

        if (facei < mesh.nInternalFaces())
        {
            const label nbri = mesh.faceNeighbour()[facei];
            const bool nbrIsCoupled = cellCoupleds[nbri];

            if (ownIsCoupled != nbrIsCoupled)
            {
                const label celli = ownIsCoupled ? owni : nbri;

                cutFacesDyn.append(facei);
                cutFacesCoupledCellDyn.append(globalCellIndex.toGlobal(celli));
            }
        }
        else
        {
            const label bFacei = facei - mesh.nInternalFaces();
            const bool nbrIsCoupled = bFaceNbrIsCoupled[bFacei];

            if (!ownIsCoupled && nbrIsCoupled)
            {
                cutFacesDyn.append(facei);
                cutFacesCoupledCellDyn.append(bFaceNbrCells[bFacei]);
            }
        }
    }

    // Transfer to fixed-size storage
    cutFaces.transfer(cutFacesDyn);
    cutFaceCoupledCells.transfer(cutFacesCoupledCellDyn);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellsToCellsExtrapolation::cellsToCellsExtrapolation()
:
    extrapolation_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellsToCellsExtrapolation::~cellsToCellsExtrapolation()
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::cellsToCellsExtrapolation>
Foam::cellsToCellsExtrapolation::New(const word& type)
{
    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(type);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName << " type "
            << type << endl << endl
            << "Valid " << typeName << " types are : " << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()();
}


// ************************************************************************* //
