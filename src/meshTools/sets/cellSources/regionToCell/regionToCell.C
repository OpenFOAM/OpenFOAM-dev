/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "regionToCell.H"
#include "regionSplit.H"
#include "emptyPolyPatch.H"
#include "cellSet.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, regionToCell, word);
    addToRunTimeSelectionTable(topoSetSource, regionToCell, istream);
}


Foam::topoSetSource::addToUsageTable Foam::regionToCell::usage_
(
    regionToCell::typeName,
    "\n    Usage: regionToCell subCellSet (pt0 .. ptn)\n\n"
    "    Select all cells in the connected region containing"
    " points (pt0..ptn).\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::regionToCell::markRegionFaces
(
    const boolList& selectedCell,
    boolList& regionFace
) const
{
    // Internal faces
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();
    forAll(faceNeighbour, facei)
    {
        if
        (
            selectedCell[faceOwner[facei]]
         != selectedCell[faceNeighbour[facei]]
        )
        {
            regionFace[facei] = true;
        }
    }

    // Swap neighbour selectedCell state
    boolList nbrSelected;
    syncTools::swapBoundaryCellList(mesh_, selectedCell, nbrSelected);

    // Boundary faces
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        const labelUList& faceCells = pp.faceCells();
        forAll(faceCells, i)
        {
            label facei = pp.start()+i;
            label bFacei = facei-mesh_.nInternalFaces();
            if (selectedCell[faceCells[i]] != nbrSelected[bFacei])
            {
                regionFace[facei] = true;
            }
        }
    }
}


Foam::boolList Foam::regionToCell::findRegions
(
    const bool verbose,
    const regionSplit& cellRegion
) const
{
    boolList keepRegion(cellRegion.nRegions(), false);

    forAll(insidePoints_, i)
    {
        // Find the region containing the insidePoint

        label celli = mesh_.findCell(insidePoints_[i]);

        label keepRegionI = -1;
        label keepProci = -1;
        if (celli != -1)
        {
            keepRegionI = cellRegion[celli];
            keepProci = Pstream::myProcNo();
        }
        reduce(keepRegionI, maxOp<label>());
        keepRegion[keepRegionI] = true;

        reduce(keepProci, maxOp<label>());

        if (keepProci == -1)
        {
            FatalErrorInFunction
                << "Did not find " << insidePoints_[i]
                << " in mesh." << " Mesh bounds are " << mesh_.bounds()
                << exit(FatalError);
        }

        if (verbose)
        {
            Info<< "    Found location " << insidePoints_[i]
                << " in cell " << celli << " on processor " << keepProci
                << " in global region " << keepRegionI
                << " out of " << cellRegion.nRegions() << " regions." << endl;
        }
    }

    return keepRegion;
}


void Foam::regionToCell::unselectOutsideRegions
(
    boolList& selectedCell
) const
{
    // Determine faces on the edge of selectedCell
    boolList blockedFace(mesh_.nFaces(), false);
    markRegionFaces(selectedCell, blockedFace);

    // Determine regions
    regionSplit cellRegion(mesh_, blockedFace);

    // Determine regions containing insidePoints_
    boolList keepRegion(findRegions(true, cellRegion));

    // Go back to bool per cell
    forAll(cellRegion, celli)
    {
        if (!keepRegion[cellRegion[celli]])
        {
            selectedCell[celli] = false;
        }
    }
}


void Foam::regionToCell::shrinkRegions
(
    boolList& selectedCell
) const
{
    // Select points on unselected cells and boundary
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    boolList boundaryPoint(mesh_.nPoints(), false);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (!pp.coupled() && !isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                const face& f = pp[i];
                forAll(f, fp)
                {
                    boundaryPoint[f[fp]] = true;
                }
            }
        }
    }

    forAll(selectedCell, celli)
    {
        if (!selectedCell[celli])
        {
            const labelList& cPoints = mesh_.cellPoints(celli);
            forAll(cPoints, i)
            {
                boundaryPoint[cPoints[i]] = true;
            }
        }
    }

    syncTools::syncPointList(mesh_, boundaryPoint, orEqOp<bool>(), false);


    // Select all cells using these points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label nChanged = 0;
    forAll(boundaryPoint, pointi)
    {
        if (boundaryPoint[pointi])
        {
            const labelList& pCells = mesh_.pointCells(pointi);
            forAll(pCells, i)
            {
                label celli = pCells[i];
                if (selectedCell[celli])
                {
                    selectedCell[celli] = false;
                    nChanged++;
                }
            }
        }
    }
    Info<< "    Eroded " << returnReduce(nChanged, sumOp<label>())
        << " cells." << endl;
}


void Foam::regionToCell::erode
(
    boolList& selectedCell
) const
{
    // Info<< "Entering shrinkRegions:" << count(selectedCell) << endl;
    // generateField("selectedCell_before", selectedCell)().write();

    // Now erode and see which regions get disconnected
    boolList shrunkSelectedCell(selectedCell);

    for (label iter = 0; iter < nErode_; iter++)
    {
        shrinkRegions(shrunkSelectedCell);
    }

    // Info<< "After shrinking:" << count(shrunkSelectedCell) << endl;
    // generateField("shrunkSelectedCell", shrunkSelectedCell)().write();



    // Determine faces on the edge of shrunkSelectedCell
    boolList blockedFace(mesh_.nFaces(), false);
    markRegionFaces(shrunkSelectedCell, blockedFace);

    // Find disconnected regions
    regionSplit cellRegion(mesh_, blockedFace);

    // Determine regions containing insidePoints
    boolList keepRegion(findRegions(true, cellRegion));


    // Extract cells in regions that are not to be kept.
    boolList removeCell(mesh_.nCells(), false);
    forAll(cellRegion, celli)
    {
        if (shrunkSelectedCell[celli] && !keepRegion[cellRegion[celli]])
        {
            removeCell[celli] = true;
        }
    }

    // Info<< "removeCell before:" << count(removeCell) << endl;
    // generateField("removeCell_before", removeCell)().write();



    // Grow removeCell
    for (label iter = 0; iter < nErode_; iter++)
    {
        // Grow selected cell in regions that are not for keeping
        boolList boundaryPoint(mesh_.nPoints(), false);
        forAll(removeCell, celli)
        {
            if (removeCell[celli])
            {
                const labelList& cPoints = mesh_.cellPoints(celli);
                forAll(cPoints, i)
                {
                    boundaryPoint[cPoints[i]] = true;
                }
            }
        }
        syncTools::syncPointList(mesh_, boundaryPoint, orEqOp<bool>(), false);

        // Select all cells using these points

        label nChanged = 0;
        forAll(boundaryPoint, pointi)
        {
            if (boundaryPoint[pointi])
            {
                const labelList& pCells = mesh_.pointCells(pointi);
                forAll(pCells, i)
                {
                    label celli = pCells[i];
                    if (!removeCell[celli])
                    {
                        removeCell[celli] = true;
                        nChanged++;
                    }
                }
            }
        }
    }

    // Info<< "removeCell after:" << count(removeCell) << endl;
    // generateField("removeCell_after", removeCell)().write();


    // Unmark removeCell
    forAll(removeCell, celli)
    {
        if (removeCell[celli])
        {
            selectedCell[celli] = false;
        }
    }
}


void Foam::regionToCell::combine(topoSet& set, const bool add) const
{
    // Note: wip. Select cells first
    boolList selectedCell(mesh_.nCells(), true);

    if (setName_.size() && setName_ != "none")
    {
        Info<< "    Loading subset " << setName_ << " to delimit search region."
            << endl;
        cellSet subSet(mesh_, setName_);

        selectedCell = false;
        forAllConstIter(cellSet, subSet, iter)
        {
            selectedCell[iter.key()] = true;
        }
    }


    unselectOutsideRegions(selectedCell);

    if (nErode_ > 0)
    {
        erode(selectedCell);
    }


    forAll(selectedCell, celli)
    {
        if (selectedCell[celli])
        {
            addOrDelete(set, celli, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionToCell::regionToCell
(
    const polyMesh& mesh,
    const word& setName,
    const pointField& insidePoints,
    const label nErode
)
:
    topoSetSource(mesh),
    setName_(setName),
    insidePoints_(insidePoints),
    nErode_(nErode)
{}


Foam::regionToCell::regionToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookupOrDefault<word>("set", "none")),
    insidePoints_
    (
        dict.found("insidePoints")
      ? dict.lookup("insidePoints")
      : dict.lookup("insidePoint")
    ),
    nErode_(dict.lookupOrDefault("nErode", 0))
{}


Foam::regionToCell::regionToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is)),
    insidePoints_(checkIs(is)),
    nErode_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionToCell::~regionToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all cells of connected region containing points "
            << insidePoints_ << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all cells of connected region containing points "
            << insidePoints_ << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
