/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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

#include "mapNearestMethod.H"
#include "nearestPatchToPatch.H"
#include "pointIndexHit.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mapNearestMethod, 0);
    addToRunTimeSelectionTable(meshToMeshMethod, mapNearestMethod, components);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::mapNearestMethod::findInitialSeeds
(
    const labelList& srcCellIDs,
    const boolList& mapFlag,
    const label startSeedI,
    label& srcSeedI,
    label& tgtSeedI
) const
{
    const vectorField& srcCcs = src_.cellCentres();

    for (label i = startSeedI; i < srcCellIDs.size(); i++)
    {
        label srcI = srcCellIDs[i];

        if (mapFlag[srcI])
        {
            const point& srcCc = srcCcs[srcI];
            pointIndexHit hit =
                tgt_.cellTree().findNearest(srcCc, great);

            if (hit.hit())
            {
                srcSeedI = srcI;
                tgtSeedI = hit.index();

                return true;
            }
            else
            {
                FatalErrorInFunction
                    << "Unable to find nearest target cell"
                    << " for source cell " << srcI
                    << " with centre " << srcCc
                    << abort(FatalError);
            }

            break;
        }
    }

    if (debug)
    {
        Pout<< "could not find starting seed" << endl;
    }

    return false;
}


Foam::scalar Foam::mapNearestMethod::calculateAddressing
(
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

    List<DynamicList<label>> srcToTgt(src_.nCells());
    List<DynamicList<label>> tgtToSrc(tgt_.nCells());

    const scalarField& srcVc = src_.cellVolumes();
    const scalarField& tgtVc = tgt_.cellVolumes();

    {
        label srcCelli = srcSeedI;
        label tgtCelli = tgtSeedI;

        do
        {
            // find nearest tgt cell
            findNearestCell(src_, tgt_, srcCelli, tgtCelli);

           // store src/tgt cell pair
            srcToTgt[srcCelli].append(tgtCelli);
            tgtToSrc[tgtCelli].append(srcCelli);

            // mark source cell srcCelli and tgtCelli as matched
            mapFlag[srcCelli] = false;

            // accumulate intersection volume
            V += srcVc[srcCelli];

            // find new source cell
            setNextNearestCells
            (
                startSeedI,
                srcCelli,
                tgtCelli,
                mapFlag,
                srcCellIDs
            );
        }
        while (srcCelli >= 0);
    }

    // for the case of multiple source cells per target cell, select the
    // nearest source cell only and discard the others
    const vectorField& srcCc = src_.cellCentres();
    const vectorField& tgtCc = tgt_.cellCentres();

    forAll(tgtToSrc, targetCelli)
    {
        if (tgtToSrc[targetCelli].size() > 1)
        {
            const vector& tgtC = tgtCc[targetCelli];

            DynamicList<label>& srcCells = tgtToSrc[targetCelli];

            label srcCelli = srcCells[0];
            scalar d = magSqr(tgtC - srcCc[srcCelli]);

            for (label i = 1; i < srcCells.size(); i++)
            {
                label srcI = srcCells[i];
                scalar dNew = magSqr(tgtC - srcCc[srcI]);
                if (dNew < d)
                {
                    d = dNew;
                    srcCelli = srcI;
                }
            }

            srcCells.clear();
            srcCells.append(srcCelli);
        }
    }

    // If there are more target cells than source cells, some target cells
    // might not yet be mapped
    forAll(tgtToSrc, tgtCelli)
    {
        if (tgtToSrc[tgtCelli].empty())
        {
            label srcCelli = findMappedSrcCell(tgtCelli, tgtToSrc);

            findNearestCell(tgt_, src_, tgtCelli, srcCelli);

            tgtToSrc[tgtCelli].append(srcCelli);
        }
    }

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


void Foam::mapNearestMethod::findNearestCell
(
    const polyMesh& mesh1,
    const polyMesh& mesh2,
    const label cell1,
    label& cell2
) const
{
    const vectorField& Cc1 = mesh1.cellCentres();
    const vectorField& Cc2 = mesh2.cellCentres();

    const vector& p1 = Cc1[cell1];

    DynamicList<label> cells2(10);
    cells2.append(cell2);

    DynamicList<label> visitedCells(10);

    scalar d = great;

    do
    {
        label c2 = cells2.remove();
        visitedCells.append(c2);

        scalar dTest = magSqr(Cc2[c2] - p1);
        if (dTest < d)
        {
            cell2 = c2;
            d = dTest;
            appendNbrCells(cell2, mesh2, visitedCells, cells2);
        }

    } while (cells2.size() > 0);
}


void Foam::mapNearestMethod::setNextNearestCells
(
    label& startSeedI,
    label& srcCelli,
    label& tgtCelli,
    boolList& mapFlag,
    const labelList& srcCellIDs
) const
{
    const labelList& srcNbr = src_.cellCells()[srcCelli];

    srcCelli = -1;
    forAll(srcNbr, i)
    {
        label celli = srcNbr[i];
        if (mapFlag[celli])
        {
            srcCelli = celli;
            return;
        }
    }

    for (label i = startSeedI; i < srcCellIDs.size(); i++)
    {
        label celli = srcCellIDs[i];
        if (mapFlag[celli])
        {
            startSeedI = i;
            break;
        }
    }

    (void)findInitialSeeds
    (
        srcCellIDs,
        mapFlag,
        startSeedI,
        srcCelli,
        tgtCelli
    );
}


Foam::label Foam::mapNearestMethod::findMappedSrcCell
(
    const label tgtCelli,
    const List<DynamicList<label>>& tgtToSrc
) const
{
    DynamicList<label> testCells(10);
    DynamicList<label> visitedCells(10);

    testCells.append(tgtCelli);

    do
    {
        // search target tgtCelli neighbours for match with source cell
        label tgtI = testCells.remove();

        if (findIndex(visitedCells, tgtI) == -1)
        {
            visitedCells.append(tgtI);

            if (tgtToSrc[tgtI].size())
            {
                return tgtToSrc[tgtI][0];
            }
            else
            {
                const labelList& nbrCells = tgt_.cellCells()[tgtI];

                forAll(nbrCells, i)
                {
                    if (findIndex(visitedCells, nbrCells[i]) == -1)
                    {
                        testCells.append(nbrCells[i]);
                    }
                }
            }
        }
    } while (testCells.size());

    // did not find any match - should not be possible to get here!
    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mapNearestMethod::mapNearestMethod
(
    const polyMesh& src,
    const polyMesh& tgt
)
:
    meshToMeshMethod(src, tgt)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::mapNearestMethod::~mapNearestMethod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::mapNearestMethod::patchToPatchMethod() const
{
    return patchToPatches::nearest::typeName;
}


Foam::scalar Foam::mapNearestMethod::calculate
(
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght,
    labelListList& tgtToSrcAddr,
    scalarListList& tgtToSrcWght
)
{
    bool ok = initialise
    (
        srcToTgtAddr,
        srcToTgtWght,
        tgtToSrcAddr,
        tgtToSrcWght
    );

    if (!ok)
    {
        return 0;
    }

    // (potentially) participating source mesh cells
    const labelList srcCellIDs(maskCells());

    // list to keep track of whether src cell can be mapped
    boolList mapFlag(src_.nCells(), false);
    UIndirectList<bool>(mapFlag, srcCellIDs) = true;

    // find initial point in tgt mesh
    label srcSeedI = -1;
    label tgtSeedI = -1;
    label startSeedI = 0;

    bool startWalk =
        findInitialSeeds
        (
            srcCellIDs,
            mapFlag,
            startSeedI,
            srcSeedI,
            tgtSeedI
        );

    if (startWalk)
    {
        const scalar V =
            calculateAddressing
            (
                srcToTgtAddr,
                srcToTgtWght,
                tgtToSrcAddr,
                tgtToSrcWght,
                srcSeedI,
                tgtSeedI,
                srcCellIDs,
                mapFlag,
                startSeedI
            );

        if (debug > 1)
        {
            writeConnectivity(src_, tgt_, srcToTgtAddr);
        }

        return V;
    }
    else
    {
        return 0;
    }
}


void Foam::mapNearestMethod::normalise
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


// ************************************************************************* //
