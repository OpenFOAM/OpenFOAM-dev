/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2024 OpenFOAM Foundation
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

#include "matchingCellsToCells.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellsToCellss
{
    defineTypeNameAndDebug(matching, 0);
    addToRunTimeSelectionTable(cellsToCells, matching, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::cellsToCellss::matching::intersect
(
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh,
    const label srcCelli,
    const label tgtCelli
) const
{
    return tgtMesh.pointInCell
    (
        srcMesh.cellCentres()[srcCelli],
        tgtCelli,
        polyMesh::FACE_PLANES
    );
}


bool Foam::cellsToCellss::matching::findInitialSeeds
(
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh,
    const labelList& srcCellIDs,
    const boolList& mapFlag,
    const label startSeedI,
    label& srcSeedI,
    label& tgtSeedI
) const
{
    const cellList& srcCells = srcMesh.cells();
    const faceList& srcFaces = srcMesh.faces();
    const pointField& srcPts = srcMesh.points();

    for (label i = startSeedI; i < srcCellIDs.size(); i++)
    {
        label srcI = srcCellIDs[i];

        if (mapFlag[srcI])
        {
            const point srcCtr(srcCells[srcI].centre(srcPts, srcFaces));
            label tgtI = tgtMesh.cellTree().findInside(srcCtr);

            if (tgtI != -1 && intersect(srcMesh, tgtMesh, srcI, tgtI))
            {
                srcSeedI = srcI;
                tgtSeedI = tgtI;

                return true;
            }
        }
    }

    if (debug)
    {
        Pout<< "could not find starting seed" << endl;
    }

    return false;
}


Foam::scalar Foam::cellsToCellss::matching::calculateAddressing
(
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh,
    labelListList& srcToTgtCellAddr,
    scalarListList& srcToTgtCellWght,
    labelListList& tgtToSrcCellAddr,
    scalarListList& tgtToSrcCellWght,
    const label srcSeedI,
    const label tgtSeedI,
    const labelList& srcCellIDs, // not used
    boolList& mapFlag,
    label& startSeedI
)
{
    scalar V = 0;

    // store a list of src cells already mapped
    labelList srcTgtSeed(srcMesh.nCells(), -1);

    List<DynamicList<label>> srcToTgt(srcMesh.nCells());
    List<DynamicList<label>> tgtToSrc(tgtMesh.nCells());

    DynamicList<label> srcSeeds(10);

    const scalarField& srcVc = srcMesh.cellVolumes();
    const scalarField& tgtVc = tgtMesh.cellVolumes();

    label srcCelli = srcSeedI;
    label tgtCelli = tgtSeedI;

    do
    {
        // store src/tgt cell pair
        srcToTgt[srcCelli].append(tgtCelli);
        tgtToSrc[tgtCelli].append(srcCelli);

        // mark source cell srcSeedI as matched
        mapFlag[srcCelli] = false;

        // accumulate intersection volume
        V += srcVc[srcCelli];

        // find new source seed cell
        appendToDirectSeeds
        (
            srcMesh,
            tgtMesh,
            mapFlag,
            srcTgtSeed,
            srcSeeds,
            srcCelli,
            tgtCelli
        );
    }
    while (srcCelli >= 0);

    // transfer addressing into persistent storage
    forAll(srcToTgtCellAddr, i)
    {
        srcToTgtCellWght[i] = scalarList(srcToTgt[i].size(), srcVc[i]);
        srcToTgtCellAddr[i].transfer(srcToTgt[i]);
    }

    forAll(tgtToSrcCellAddr, i)
    {
        tgtToSrcCellWght[i] = scalarList(tgtToSrc[i].size(), tgtVc[i]);
        tgtToSrcCellAddr[i].transfer(tgtToSrc[i]);
    }

    return V;
}


void Foam::cellsToCellss::matching::appendToDirectSeeds
(
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh,
    boolList& mapFlag,
    labelList& srcTgtSeed,
    DynamicList<label>& srcSeeds,
    label& srcSeedI,
    label& tgtSeedI
) const
{
    const labelList& srcNbr = srcMesh.cellCells()[srcSeedI];
    const labelList& tgtNbr = tgtMesh.cellCells()[tgtSeedI];

    forAll(srcNbr, i)
    {
        label srcI = srcNbr[i];

        if (mapFlag[srcI] && (srcTgtSeed[srcI] == -1))
        {
            // source cell srcI not yet mapped

            // identify if target cell exists for source cell srcI
            bool found = false;
            forAll(tgtNbr, j)
            {
                label tgtI = tgtNbr[j];

                if (intersect(srcMesh, tgtMesh, srcI, tgtI))
                {
                    // new match - append to lists
                    found = true;

                    srcTgtSeed[srcI] = tgtI;
                    srcSeeds.append(srcI);

                    break;
                }
            }

            if (!found)
            {
                // no match available for source cell srcI
                mapFlag[srcI] = false;
            }
        }
    }

    if (srcSeeds.size())
    {
        srcSeedI = srcSeeds.remove();
        tgtSeedI = srcTgtSeed[srcSeedI];
    }
    else
    {
        srcSeedI = -1;
        tgtSeedI = -1;
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::cellsToCellss::matching::calculate
(
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh
)
{
    initialise(srcMesh, tgtMesh);

    // Determine (potentially) participating source mesh cells
    const labelList srcCellIDs(maskCells(srcMesh, tgtMesh));

    // Initialise list to keep track of whether src cell can be mapped
    boolList mapFlag(srcMesh.nCells(), false);
    UIndirectList<bool>(mapFlag, srcCellIDs) = true;

    // Find initial point in tgt mesh
    label srcSeedI = -1;
    label tgtSeedI = -1;
    label startSeedI = 0;
    bool startWalk =
        findInitialSeeds
        (
            srcMesh,
            tgtMesh,
            srcCellIDs,
            mapFlag,
            startSeedI,
            srcSeedI,
            tgtSeedI
        );

    if (startWalk)
    {
        return
            calculateAddressing
            (
                srcMesh,
                tgtMesh,
                srcLocalTgtCells_,
                srcWeights_,
                tgtLocalSrcCells_,
                tgtWeights_,
                srcSeedI,
                tgtSeedI,
                srcCellIDs,
                mapFlag,
                startSeedI
            );
    }
    else
    {
        return 0;
    }
}


void Foam::cellsToCellss::matching::normalise
(
    const polyMesh& srcMesh,
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght
) const
{
    forAll(srcToTgtWght, srcCelli)
    {
        if (srcToTgtWght[srcCelli].size() > 1)
        {
            srcToTgtAddr[srcCelli].resize(1);
            srcToTgtWght[srcCelli].resize(1);
            srcToTgtWght[srcCelli][0] = 1;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellsToCellss::matching::matching()
:
    cellsToCells()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellsToCellss::matching::~matching()
{}


// ************************************************************************* //
