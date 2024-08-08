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

#include "scotch.H"
#include "addToRunTimeSelectionTable.H"
#include "floatScalar.H"
#include "Time.H"
#include "OFstream.H"
#include "globalIndex.H"
#include "SubField.H"

extern "C"
{
#include "scotch.h"
}


// Hack: scotch generates floating point errors so need to switch of error
//       trapping!
#ifdef __GLIBC__
    #ifndef _GNU_SOURCE
        #define _GNU_SOURCE
    #endif
    #include <fenv.h>
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionMethods
{
    defineTypeNameAndDebug(scotch, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        scotch,
        decomposer
    );

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        scotch,
        distributor
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::decompositionMethods::scotch::check
(
    const int retVal,
    const char* str
)
{
    if (retVal)
    {
        FatalErrorInFunction
            << "Call to scotch routine " << str << " failed."
            << exit(FatalError);
    }
}


Foam::label Foam::decompositionMethods::scotch::decompose
(
    const fileName& meshPath,
    const labelList& adjncy,
    const labelList& xadj,
    const scalarField& cellWeights,
    labelList& decomp
)
{
    if (!Pstream::parRun())
    {
        decomposeOneProc
        (
            meshPath,
            adjncy,
            xadj,
            cellWeights,
            decomp
        );
    }
    else
    {
        if (debug)
        {
            Info<< "scotch : running in parallel."
                << " Decomposing all of graph on master processor." << endl;
        }
        globalIndex globalCells(xadj.size()-1);
        label nTotalConnections = returnReduce(adjncy.size(), sumOp<label>());

        // Send all to master. Use scheduled to save some storage.
        if (Pstream::master())
        {
            Field<label> allAdjncy(nTotalConnections);
            Field<label> allXadj(globalCells.size()+1);
            scalarField allWeights(globalCells.size());

            // Insert my own
            label nTotalCells = 0;
            forAll(cellWeights, celli)
            {
                allXadj[nTotalCells] = xadj[celli];
                allWeights[nTotalCells++] = cellWeights[celli];
            }
            nTotalConnections = 0;
            forAll(adjncy, i)
            {
                allAdjncy[nTotalConnections++] = adjncy[i];
            }

            for (int slave=1; slave<Pstream::nProcs(); slave++)
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                Field<label> nbrAdjncy(fromSlave);
                Field<label> nbrXadj(fromSlave);
                scalarField nbrWeights(fromSlave);

                // Append.
                // label procStart = nTotalCells;
                forAll(nbrXadj, celli)
                {
                    allXadj[nTotalCells] = nTotalConnections+nbrXadj[celli];
                    allWeights[nTotalCells++] = nbrWeights[celli];
                }
                // No need to renumber xadj since already global.
                forAll(nbrAdjncy, i)
                {
                    allAdjncy[nTotalConnections++] = nbrAdjncy[i];
                }
            }
            allXadj[nTotalCells] = nTotalConnections;


            Field<label> allFinalDecomp;
            decomposeOneProc
            (
                meshPath,
                allAdjncy,
                allXadj,
                allWeights,
                allFinalDecomp
            );


            // Send allFinalDecomp back
            for (int slave=1; slave<Pstream::nProcs(); slave++)
            {
                OPstream toSlave(Pstream::commsTypes::scheduled, slave);
                toSlave << SubField<label>
                (
                    allFinalDecomp,
                    globalCells.localSize(slave),
                    globalCells.offset(slave)
                );
            }
            // Get my own part (always first)
            decomp = SubField<label>
            (
                allFinalDecomp,
                globalCells.localSize()
            );
        }
        else
        {
            // Send my part of the graph (already in global numbering)
            {
                OPstream toMaster
                (
                    Pstream::commsTypes::scheduled,
                    Pstream::masterNo()
                );
                toMaster<< adjncy << SubField<label>(xadj, xadj.size()-1)
                    << cellWeights;
            }

            // Receive back decomposition
            IPstream fromMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );
            fromMaster >> decomp;
        }
    }
    return 0;
}


Foam::label Foam::decompositionMethods::scotch::decomposeOneProc
(
    const fileName& meshPath,
    const labelList& adjncy,
    const labelList& xadj,
    const scalarField& cellWeights,

    labelList& decomp
)
{
    // Dump graph
    if (!methodDict_.empty())
    {
        if (methodDict_.lookupOrDefault("writeGraph", false))
        {
            OFstream str(meshPath + ".grf");

            Info<< "Dumping Scotch graph file to " << str.name() << endl
                << "Use this in combination with gpart." << endl;

            label version = 0;
            str << version << nl;
            // Number of vertices
            str << xadj.size()-1 << ' ' << adjncy.size() << nl;
            // Numbering starts from 0
            label baseval = 0;
            // Has weights?
            label hasEdgeWeights = 0;
            label hasVertexWeights = 0;
            label numericflag = 10*hasEdgeWeights+hasVertexWeights;
            str << baseval << ' ' << numericflag << nl;
            for (label celli = 0; celli < xadj.size()-1; celli++)
            {
                label start = xadj[celli];
                label end = xadj[celli+1];
                str << end-start;

                for (label i = start; i < end; i++)
                {
                    str << ' ' << adjncy[i];
                }
                str << nl;
            }
        }
    }

    // Reset the seed of the pseudo-random generator used by the graph
    // partitioning routines of the libScotch library. Two consecutive calls to
    // the same libScotch partitioning routines, and separated by a call to
    // SCOTCH randomReset, will always yield the same results, as if the
    // equivalent standalone Scotch programs were used twice, independently,
    SCOTCH_randomReset();

    // Strategy
    // ~~~~~~~~

    // Default.
    SCOTCH_Strat stradat;
    check(SCOTCH_stratInit(&stradat), "SCOTCH_stratInit");

    if (!methodDict_.empty())
    {
        string strategy;
        if (methodDict_.readIfPresent("strategy", strategy))
        {
            if (debug)
            {
                Info<< "scotch : Using strategy " << strategy << endl;
            }
            SCOTCH_stratGraphMap(&stradat, strategy.c_str());
        }
    }


    // Graph
    // ~~~~~

    labelList velotab;

    // Check for externally provided cellweights and if so initialise weights
    if (!cellWeights.empty())
    {
        if (cellWeights.size() != xadj.size()-1)
        {
            FatalErrorInFunction
                << "Number of cell weights " << cellWeights.size()
                << " does not equal number of cells " << xadj.size()-1
                << exit(FatalError);
        }

        label nWeights = 1;
        velotab = scaleWeights(cellWeights, nWeights, false);
    }


    SCOTCH_Graph grafdat;
    check(SCOTCH_graphInit(&grafdat), "SCOTCH_graphInit");
    check
    (
        SCOTCH_graphBuild
        (
            &grafdat,
            0,                      // baseval, c-style numbering
            xadj.size()-1,          // vertnbr, nCells
            xadj.begin(),           // verttab, start index per cell into adjncy
            &xadj[1],               // vendtab, end index  ,,
            velotab.begin(),        // velotab, vertex weights
            nullptr,                // vlbltab
            adjncy.size(),          // edgenbr, number of arcs
            adjncy.begin(),         // edgetab
            nullptr                 // edlotab, edge weights
        ),
        "SCOTCH_graphBuild"
    );
    check(SCOTCH_graphCheck(&grafdat), "SCOTCH_graphCheck");


    // Architecture
    // ~~~~~~~~~~~~
    // (fully connected network topology since using switch)

    SCOTCH_Arch archdat;
    check(SCOTCH_archInit(&archdat), "SCOTCH_archInit");

    labelList processorWeights;
    if (!methodDict_.empty())
    {
        methodDict_.readIfPresent("processorWeights", processorWeights);
    }
    if (processorWeights.size())
    {
        if (debug)
        {
            Info<< "scotch : Using processor weights " << processorWeights
                << endl;
        }
        check
        (
            SCOTCH_archCmpltw(&archdat, nProcessors_, processorWeights.begin()),
            "SCOTCH_archCmpltw"
        );
    }
    else
    {
        check
        (
            SCOTCH_archCmplt(&archdat, nProcessors_),
            "SCOTCH_archCmplt"
        );
    }

    // Hack:switch off fpu error trapping
    #ifdef FE_NOMASK_ENV
    int oldExcepts = fedisableexcept
    (
        FE_DIVBYZERO
      | FE_INVALID
      | FE_OVERFLOW
    );
    #endif

    decomp.setSize(xadj.size()-1);
    decomp = 0;
    check
    (
        SCOTCH_graphMap
        (
            &grafdat,
            &archdat,
            &stradat,           // const SCOTCH_Strat *
            decomp.begin() // parttab
        ),
        "SCOTCH_graphMap"
    );

    #ifdef FE_NOMASK_ENV
    feenableexcept(oldExcepts);
    #endif

    // Release storage for graph
    SCOTCH_graphExit(&grafdat);

    // Release storage for strategy
    SCOTCH_stratExit(&stradat);

    // Release storage for network topology
    SCOTCH_archExit(&archdat);

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionMethods::scotch::scotch
(
    const dictionary& decompositionDict,
    const dictionary& methodDict
)
:
    decompositionMethod(decompositionDict),
    methodDict_(methodDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::decompositionMethods::scotch::decompose
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

    // Calculate local or global (if Pstream::parRun()) connectivity
    CompactListList<label> cellCells;
    calcCellCells
    (
        mesh,
        identityMap(mesh.nCells()),
        mesh.nCells(),
        true,
        cellCells
    );

    // Decompose using default weights
    labelList decomp;
    decompose
    (
        mesh.time().path()/mesh.name(),
        cellCells.m(),
        cellCells.offsets(),
        pointWeights,
        decomp
    );

    return decomp;
}


Foam::labelList Foam::decompositionMethods::scotch::decompose
(
    const polyMesh& mesh,
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& pointWeights
)
{
    if (agglom.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Size of cell-to-coarse map " << agglom.size()
            << " differs from number of cells in mesh " << mesh.nCells()
            << exit(FatalError);
    }

    checkWeights(agglomPoints, pointWeights);

    // Calculate local or global (if Pstream::parRun()) connectivity
    CompactListList<label> cellCells;
    calcCellCells
    (
        mesh,
        agglom,
        agglomPoints.size(),
        true,
        cellCells
    );

    // Decompose using weights
    labelList decomp;
    decompose
    (
        mesh.time().path()/mesh.name(),
        cellCells.m(),
        cellCells.offsets(),
        pointWeights,
        decomp
    );

    // Rework back into decomposition for original mesh_
    labelList fineDistribution(agglom.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = decomp[agglom[i]];
    }

    return fineDistribution;
}


Foam::labelList Foam::decompositionMethods::scotch::decompose
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

    // Decompose using weights
    labelList decomp;
    decompose
    (
        "scotch",
        cellCells.m(),
        cellCells.offsets(),
        cellWeights,
        decomp
    );

    return decomp;
}


// ************************************************************************* //
