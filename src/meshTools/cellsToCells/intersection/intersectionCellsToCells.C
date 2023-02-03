/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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

#include "intersectionCellsToCells.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "tetOverlapVolume.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellsToCellss
{
    defineTypeNameAndDebug(intersection, 0);
    addToRunTimeSelectionTable(cellsToCells, intersection, word);
}
}


const Foam::scalar Foam::cellsToCellss::intersection::tolerance_ = 1e-6;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::cellsToCellss::intersection::intersect
(
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh,
    const label srcCelli,
    const label tgtCelli
) const
{
    return
        tetOverlapVolume().cellCellOverlapMinDecomp
        (
            srcMesh,
            srcCelli,
            tgtMesh,
            tgtCelli,
            treeBoundBox(tgtMesh.points(), tgtMesh.cellPoints()[tgtCelli]),
            tolerance_*srcMesh.cellVolumes()[srcCelli]
        );
}


Foam::scalar Foam::cellsToCellss::intersection::interVol
(
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh,
    const label srcCelli,
    const label tgtCelli
) const
{
    return
        tetOverlapVolume().cellCellOverlapVolumeMinDecomp
        (
            srcMesh,
            srcCelli,
            tgtMesh,
            tgtCelli,
            treeBoundBox(tgtMesh.points(), tgtMesh.cellPoints()[tgtCelli])
        );
}


bool Foam::cellsToCellss::intersection::findInitialSeeds
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
        const label srcI = srcCellIDs[i];

        if (mapFlag[srcI])
        {
            const labelList tgtIDs
            (
                tgtMesh.cellTree().findBox
                (
                    treeBoundBox(srcCells[srcI].bb(srcPts, srcFaces))
                )
            );

            forAll(tgtIDs, j)
            {
                const label tgtI = tgtIDs[j];

                if (intersect(srcMesh, tgtMesh, srcI, tgtI))
                {
                    srcSeedI = srcI;
                    tgtSeedI = tgtI;

                    return true;
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "could not find starting seed" << endl;
    }

    return false;
}


Foam::scalar Foam::cellsToCellss::intersection::calculateAddressing
(
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh,
    labelListList& srcToTgtCellAddr,
    scalarListList& srcToTgtCellWght,
    labelListList& tgtToSrcCellAddr,
    scalarListList& tgtToSrcCellWght,
    const label srcSeedI,
    const label tgtSeedI,
    const labelList& srcCellIDs,
    boolList& mapFlag,
    label& startSeedI
)
{
    scalar V = 0;

    label srcCelli = srcSeedI;
    label tgtCelli = tgtSeedI;

    List<DynamicList<label>> srcToTgtAddr(srcMesh.nCells());
    List<DynamicList<scalar>> srcToTgtWght(srcMesh.nCells());

    List<DynamicList<label>> tgtToSrcAddr(tgtMesh.nCells());
    List<DynamicList<scalar>> tgtToSrcWght(tgtMesh.nCells());

    // list of tgt cell neighbour cells
    DynamicList<label> nbrTgtCells(10);

    // list of tgt cells currently visited for srcCelli to avoid multiple hits
    DynamicList<label> visitedTgtCells(10);

    // list to keep track of tgt cells used to seed src cells
    labelList seedCells(srcMesh.nCells(), -1);
    seedCells[srcCelli] = tgtCelli;

    const scalarField& srcVol = srcMesh.cellVolumes();

    do
    {
        nbrTgtCells.clear();
        visitedTgtCells.clear();

        // append initial target cell and neighbours
        nbrTgtCells.append(tgtCelli);
        appendNbrCells(tgtCelli, tgtMesh, visitedTgtCells, nbrTgtCells);

        do
        {
            tgtCelli = nbrTgtCells.remove();
            visitedTgtCells.append(tgtCelli);

            scalar vol = interVol(srcMesh, tgtMesh, srcCelli, tgtCelli);

            // accumulate addressing and weights for valid intersection
            if (vol/srcVol[srcCelli] > tolerance_)
            {
                // store src/tgt cell pair
                srcToTgtAddr[srcCelli].append(tgtCelli);
                srcToTgtWght[srcCelli].append(vol);

                tgtToSrcAddr[tgtCelli].append(srcCelli);
                tgtToSrcWght[tgtCelli].append(vol);

                appendNbrCells(tgtCelli, tgtMesh, visitedTgtCells, nbrTgtCells);

                // accumulate intersection volume
                V += vol;
            }
        }
        while (!nbrTgtCells.empty());

        mapFlag[srcCelli] = false;

        // find new source seed cell
        setNextCells
        (
            srcMesh,
            tgtMesh,
            startSeedI,
            srcCelli,
            tgtCelli,
            srcCellIDs,
            mapFlag,
            visitedTgtCells,
            seedCells
        );
    }
    while (srcCelli != -1);

    // transfer addressing into persistent storage
    forAll(srcToTgtCellAddr, i)
    {
        srcToTgtCellAddr[i].transfer(srcToTgtAddr[i]);
        srcToTgtCellWght[i].transfer(srcToTgtWght[i]);
    }

    forAll(tgtToSrcCellAddr, i)
    {
        tgtToSrcCellAddr[i].transfer(tgtToSrcAddr[i]);
        tgtToSrcCellWght[i].transfer(tgtToSrcWght[i]);
    }

    return V;
}


void Foam::cellsToCellss::intersection::setNextCells
(
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh,
    label& startSeedI,
    label& srcCelli,
    label& tgtCelli,
    const labelList& srcCellIDs,
    const boolList& mapFlag,
    const DynamicList<label>& visitedCells,
    labelList& seedCells
) const
{
    const labelList& srcNbrCells = srcMesh.cellCells()[srcCelli];

    // set possible seeds for later use by querying all src cell neighbours
    // with all visited target cells
    bool valuesSet = false;
    forAll(srcNbrCells, i)
    {
        label cellS = srcNbrCells[i];

        if (mapFlag[cellS] && seedCells[cellS] == -1)
        {
            forAll(visitedCells, j)
            {
                label cellT = visitedCells[j];

                if (intersect(srcMesh, tgtMesh, cellS, cellT))
                {
                    seedCells[cellS] = cellT;

                    if (!valuesSet)
                    {
                        srcCelli = cellS;
                        tgtCelli = cellT;
                        valuesSet = true;
                    }
                }
            }
        }
    }

    // set next src and tgt cells if not set above
    if (valuesSet)
    {
        return;
    }
    else
    {
        // try to use existing seed
        bool foundNextSeed = false;
        for (label i = startSeedI; i < srcCellIDs.size(); i++)
        {
            label cellS = srcCellIDs[i];

            if (mapFlag[cellS])
            {
                if (!foundNextSeed)
                {
                    startSeedI = i;
                    foundNextSeed = true;
                }

                if (seedCells[cellS] != -1)
                {
                    srcCelli = cellS;
                    tgtCelli = seedCells[cellS];

                    return;
                }
            }
        }

        // perform new search to find match
        if (debug)
        {
            Pout<< "Advancing front stalled: searching for new "
                << "target cell" << endl;
        }

        bool restart =
            findInitialSeeds
            (
                srcMesh,
                tgtMesh,
                srcCellIDs,
                mapFlag,
                startSeedI,
                srcCelli,
                tgtCelli
            );

        if (restart)
        {
            // successfully found new starting seed-pair
            return;
        }
    }

    // if we have got to here, there are no more src/tgt cell intersections
    srcCelli = -1;
    tgtCelli = -1;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::cellsToCellss::intersection::calculate
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


void Foam::cellsToCellss::intersection::normalise
(
    const polyMesh& srcMesh,
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght
) const
{
    tetOverlapVolume overlapEngine;

    forAll(srcToTgtWght, srcCelli)
    {
        const scalar v = overlapEngine.cellVolumeMinDecomp(srcMesh, srcCelli);

        forAll(srcToTgtWght[srcCelli], i)
        {
            srcToTgtWght[srcCelli][i] /= v;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellsToCellss::intersection::intersection()
:
    cellsToCells()
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::cellsToCellss::intersection::~intersection()
{}


// ************************************************************************* //
