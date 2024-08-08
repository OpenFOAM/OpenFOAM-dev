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

#include "metis.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

extern "C"
{
    #define OMPI_SKIP_MPICXX
    #include "metis.h"
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionMethods
{
    defineTypeNameAndDebug(metis, 0);
    addToRunTimeSelectionTable(decompositionMethod, metis, decomposer);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::decompositionMethods::metis::decompose
(
    const labelList& adjncy,
    const labelList& xadj,
    const scalarField& cellWeights,
    labelList& decomp
)
{
    // Method of decomposition
    // recursive: multi-level recursive bisection (default)
    // kWay: multi-level k-way
    word method("recursive");

    label numCells = xadj.size()-1;

    // Decomposition options
    labelList options(METIS_NOPTIONS);
    METIS_SetDefaultOptions(options.begin());

    // Processor weights initialised with no size, only used if specified in
    // a file
    Field<real_t> processorWeights;

    if (cellWeights.size() > 0 && cellWeights.size() != numCells)
    {
        FatalErrorInFunction
            << "Number of cell weights " << cellWeights.size()
            << " does not equal number of cells " << numCells
            << exit(FatalError);
    }

    label nWeights = 1;

    // Cell weights (so on the vertices of the dual)
    labelList intWeights(scaleWeights(cellWeights, nWeights, false));

    // Face weights (so on the edges of the dual)
    labelList faceWeights;


    // Check for user supplied weights and decomp options
    if (!methodDict_.empty())
    {
        if (methodDict_.readIfPresent("method", method))
        {
            if (method != "recursive" && method != "kWay")
            {
                FatalIOErrorInFunction(methodDict_)
                    << "Method " << method
                    << " should be 'recursive' or 'kWay'"
                    << exit(FatalIOError);
            }

            Info<< "metis : Using Metis method     " << method
                << nl << endl;
        }

        if (methodDict_.readIfPresent("options", options))
        {
            if (options.size() != METIS_NOPTIONS)
            {
                FatalIOErrorInFunction(methodDict_)
                    << "Number of options should be " << METIS_NOPTIONS
                    << " found " << options
                    << exit(FatalIOError);
            }

            Info<< "Using Metis options     " << options << nl << endl;
        }

        if (methodDict_.readIfPresent("processorWeights", processorWeights))
        {
            processorWeights /= sum(processorWeights);

            if (processorWeights.size() != nProcessors_)
            {
                FatalIOErrorInFunction(methodDict_)
                    << "Number of processor weights "
                    << processorWeights.size()
                    << " does not equal number of domains " << nProcessors_
                    << exit(FatalIOError);
            }
        }
    }

    label ncon = 1;
    label nProcs = nProcessors_;

    // Output: cell -> processor addressing
    decomp.setSize(numCells);

    // Output: number of cut edges
    label edgeCut = 0;

    if (method == "recursive")
    {
        METIS_PartGraphRecursive
        (
            &numCells,                  // num vertices in graph
            &ncon,                      // num balancing constraints
            const_cast<labelList&>(xadj).begin(),   // indexing into adjncy
            const_cast<labelList&>(adjncy).begin(), // neighbour info
            intWeights.begin(),         // vertexweights
            nullptr,                    // vsize: total communication vol
            faceWeights.begin(),        // edgeweights
            &nProcs,                    // nParts
            processorWeights.begin(),   // tpwgts
            nullptr,                    // ubvec: processor imbalance (default)
            options.begin(),
            &edgeCut,
            decomp.begin()
        );
    }
    else
    {
        METIS_PartGraphKway
        (
            &numCells,                  // num vertices in graph
            &ncon,                      // num balancing constraints
            const_cast<labelList&>(xadj).begin(),   // indexing into adjncy
            const_cast<labelList&>(adjncy).begin(), // neighbour info
            intWeights.begin(),         // vertexweights
            nullptr,                    // vsize: total communication vol
            faceWeights.begin(),        // edgeweights
            &nProcs,                    // nParts
            processorWeights.begin(),   // tpwgts
            nullptr,                    // ubvec: processor imbalance (default)
            options.begin(),
            &edgeCut,
            decomp.begin()
        );
    }

    return edgeCut;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionMethods::metis::metis
(
    const dictionary& decompositionDict,
    const dictionary& methodDict
)
:
    decompositionMethod(decompositionDict),
    methodDict_(methodDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::decompositionMethods::metis::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
)
{
    if (points.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Can use this decomposition method only for the whole mesh"
            << endl
            << "and supply one coordinate (cellCentre) for every cell." << endl
            << "The number of coordinates " << points.size() << endl
            << "The number of cells in the mesh " << mesh.nCells()
            << exit(FatalError);
    }

    checkWeights(points, pointWeights);

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
    CompactListList<label> cellCells;
    calcCellCells
    (
        mesh,
        identityMap(mesh.nCells()),
        mesh.nCells(),
        false,
        cellCells
    );

    // Decompose using default weights
    labelList decomp;
    decompose
    (
        cellCells.m(),
        cellCells.offsets(),
        pointWeights,
        decomp
    );

    return decomp;
}


Foam::labelList Foam::decompositionMethods::metis::decompose
(
    const polyMesh& mesh,
    const labelList& cellToRegion,
    const pointField& regionPoints,
    const scalarField& regionWeights
)
{
    if (cellToRegion.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Size of cell-to-coarse map " << cellToRegion.size()
            << " differs from number of cells in mesh " << mesh.nCells()
            << exit(FatalError);
    }

    checkWeights(regionPoints, regionWeights);

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    CompactListList<label> cellCells;
    calcCellCells(mesh, cellToRegion, regionPoints.size(), false, cellCells);

    // Decompose using default weights
    labelList decomp;
    decompose(cellCells.m(), cellCells.offsets(), regionWeights, decomp);


    // Rework back into decomposition for original mesh
    labelList fineDistribution(cellToRegion.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = decomp[cellToRegion[i]];
    }

    return fineDistribution;
}


Foam::labelList Foam::decompositionMethods::metis::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cellWeights
)
{
    if (cellCentres.size() != globalCellCells.size())
    {
        FatalErrorInFunction
            << "Inconsistent number of cells (" << globalCellCells.size()
            << ") and number of cell centres (" << cellCentres.size()
            << ")." << exit(FatalError);
    }

    checkWeights(cellCentres, cellWeights);

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
    CompactListList<label> cellCells(globalCellCells);

    // Decompose using default weights
    labelList decomp;
    decompose(cellCells.m(), cellCells.offsets(), cellWeights, decomp);

    return decomp;
}


// ************************************************************************* //
