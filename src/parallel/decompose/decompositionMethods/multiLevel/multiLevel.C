/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "multiLevel.H"
#include "IFstream.H"
#include "globalIndex.H"
#include "distributionMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionMethods
{
    defineTypeNameAndDebug(multiLevel, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        multiLevel,
        decomposer
    );

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        multiLevel,
        distributor
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::decompositionMethods::multiLevel::subsetGlobalCellCells
(
    const label nDomains,
    const label domainI,
    const labelList& dist,

    const labelListList& cellCells,
    const labelList& set,
    labelListList& subCellCells,
    labelList& cutConnections
) const
{
    // Determine new index for cells by inverting subset
    labelList oldToNew(invert(cellCells.size(), set));

    globalIndex globalCells(cellCells.size());

    // Subset locally the elements for which I have data
    subCellCells = UIndirectList<labelList>(cellCells, set);

    // Get new indices for neighbouring processors
    List<Map<label>> compactMap;
    distributionMap map(globalCells, subCellCells, compactMap);
    map.distribute(oldToNew);
    labelList allDist(dist);
    map.distribute(allDist);

    // Now we have:
    // oldToNew : the locally-compact numbering of all our cellCells. -1 if
    //            cellCell is not in set.
    // allDist  : destination domain for all our cellCells
    // subCellCells : indexes into oldToNew and allDist

    // Globally compact numbering for cells in set.
    globalIndex globalSubCells(set.size());

    // Now subCellCells contains indices into oldToNew which are the
    // new locations of the neighbouring cells.

    cutConnections.setSize(nDomains);
    cutConnections = 0;

    forAll(subCellCells, subCelli)
    {
        labelList& cCells = subCellCells[subCelli];

        // Keep the connections to valid mapped cells
        label newI = 0;
        forAll(cCells, i)
        {
            // Get locally-compact cell index of neighbouring cell
            label nbrCelli = oldToNew[cCells[i]];
            if (nbrCelli == -1)
            {
                cutConnections[allDist[cCells[i]]]++;
            }
            else
            {
                // Reconvert local cell index into global one

                // Get original neighbour
                label celli = set[subCelli];
                label oldNbrCelli = cellCells[celli][i];
                // Get processor from original neighbour
                label proci = globalCells.whichProcID(oldNbrCelli);
                // Convert into global compact numbering
                cCells[newI++] = globalSubCells.toGlobal(proci, nbrCelli);
            }
        }
        cCells.setSize(newI);
    }
}


void Foam::decompositionMethods::multiLevel::decompose
(
    const labelListList& pointPoints,
    const pointField& points,
    const scalarField& pointWeights,
    const labelList& pointMap,      // map back to original points
    const label levelI,
    labelField& decomp
)
{
    labelList dist
    (
        methods_[levelI].decompose
        (
            pointPoints,
            points,
            pointWeights
        )
    );

    forAll(pointMap, i)
    {
        label orig = pointMap[i];
        decomp[orig] += dist[i];
    }

    if (levelI != methods_.size()-1)
    {
        // Recurse

        // Determine points per domain
        label n = methods_[levelI].nDomains();
        labelListList domainToPoints(invertOneToMany(n, dist));

        // 'Make space' for new levels of decomposition
        decomp *= methods_[levelI+1].nDomains();

        // Extract processor+local index from point-point addressing
        if (debug && Pstream::master())
        {
            Pout<< "Decomposition at level " << levelI << " :" << endl;
        }

        for (label domainI = 0; domainI < n; domainI++)
        {
            // Extract elements for current domain
            const labelList domainPoints(findIndices(dist, domainI));

            // Subset point-wise data.
            pointField subPoints(points, domainPoints);
            scalarField subWeights(pointWeights, domainPoints);
            labelList subPointMap(UIndirectList<label>(pointMap, domainPoints));
            // Subset point-point addressing (adapt global numbering)
            labelListList subPointPoints;
            labelList nOutsideConnections;
            subsetGlobalCellCells
            (
                n,
                domainI,
                dist,

                pointPoints,
                domainPoints,

                subPointPoints,
                nOutsideConnections
            );

            label nPoints = returnReduce(domainPoints.size(), plusOp<label>());
            Pstream::listCombineGather(nOutsideConnections, plusEqOp<label>());
            Pstream::listCombineScatter(nOutsideConnections);
            label nPatches = 0;
            label nFaces = 0;
            forAll(nOutsideConnections, i)
            {
                if (nOutsideConnections[i] > 0)
                {
                    nPatches++;
                    nFaces += nOutsideConnections[i];
                }
            }

            string oldPrefix;
            if (debug && Pstream::master())
            {
                Pout<< "    Domain " << domainI << nl
                    << "        Number of cells = " << nPoints << nl
                    << "        Number of inter-domain patches = " << nPatches
                    << nl
                    << "        Number of inter-domain faces = " << nFaces << nl
                    << endl;
                oldPrefix = Pout.prefix();
                Pout.prefix() = "  " + oldPrefix;
            }

            decompose
            (
                subPointPoints,
                subPoints,
                subWeights,
                subPointMap,
                levelI+1,

                decomp
            );
            if (debug && Pstream::master())
            {
                Pout.prefix() = oldPrefix;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionMethods::multiLevel::multiLevel
(
    const dictionary& decompositionDict,
    const dictionary& methodDict
)
:
    decompositionMethod(decompositionDict)
{
    methods_.setSize(methodDict.size());
    label i = 0;
    forAllConstIter(dictionary, methodDict, iter)
    {
        methods_.set(i++, decompositionMethod::NewDecomposer(iter().dict()));
    }

    label n = 1;
    Info<< "decompositionMethod " << type() << " :" << endl;
    forAll(methods_, i)
    {
        Info<< "    level " << i << " decomposing with " << methods_[i].type()
            << " into " << methods_[i].nDomains() << " subdomains." << endl;

        n *= methods_[i].nDomains();
    }

    if (n != nDomains())
    {
        FatalErrorInFunction
            << "Top level decomposition specifies " << nDomains()
            << " domains which is not equal to the product of"
            << " all sub domains " << n
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::decompositionMethods::multiLevel::decompose
(
    const polyMesh& mesh,
    const pointField& cellCentres,
    const scalarField& cellWeights
)
{
    CompactListList<label> cellCells;
    calcCellCells
    (
        mesh,
        identityMap(cellCentres.size()),
        cellCentres.size(),
        true,
        cellCells
    );

    labelField decomp(cellCentres.size(), 0);
    labelList cellMap(identityMap(cellCentres.size()));

    decompose
    (
        cellCells.list(),
        cellCentres,
        cellWeights,
        cellMap,      // map back to original cells
        0,

        decomp
    );

    return decomp;
}


Foam::labelList Foam::decompositionMethods::multiLevel::decompose
(
    const labelListList& globalPointPoints,
    const pointField& points,
    const scalarField& pointWeights
)
{
    labelField decomp(points.size(), 0);
    labelList pointMap(identityMap(points.size()));

    decompose
    (
        globalPointPoints,
        points,
        pointWeights,
        pointMap,       // map back to original points
        0,

        decomp
    );

    return decomp;
}


// ************************************************************************* //
