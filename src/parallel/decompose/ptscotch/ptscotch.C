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

#include "ptscotch.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "OFstream.H"
#include "globalIndex.H"
#include "SubField.H"
#include "PstreamGlobals.H"

extern "C"
{
    #include <stdio.h>
    #include <mpi.h>
    #include "ptscotch.h"
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
    defineTypeNameAndDebug(ptscotch, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        ptscotch,
        distributor
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::decompositionMethods::ptscotch::check
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


Foam::label Foam::decompositionMethods::ptscotch::decompose
(
    const fileName& meshPath,
    const labelList& adjncy,
    const labelList& xadj,
    const scalarField& cellWeights,
    labelList& decomp
) const
{
    labelList dummyAdjncy(1);
    labelList dummyXadj(1, label(0));

    return decompose
    (
        meshPath,
        adjncy.size(),
        (adjncy.size() ? adjncy.begin() : dummyAdjncy.begin()),
        xadj.size(),
        (xadj.size() ? xadj.begin() : dummyXadj.begin()),
        cellWeights,
        decomp
    );
}


Foam::label Foam::decompositionMethods::ptscotch::decompose
(
    const fileName& meshPath,
    const label adjncySize,
    const label adjncy[],
    const label xadjSize,
    const label xadj[],
    const scalarField& cellWeights,
    labelList& decomp
) const
{
    if (debug)
    {
        Pout<< "ptscotch : entering with xadj:" << xadjSize << endl;
    }

    // Dump graph
    if (decompositionDict_.found("scotchCoeffs"))
    {
        const dictionary& scotchCoeffs =
            decompositionDict_.subDict("scotchCoeffs");

        if (scotchCoeffs.lookupOrDefault("writeGraph", false))
        {
            OFstream str
            (
               meshPath + "_" + Foam::name(Pstream::myProcNo()) + ".dgr"
            );

            Pout<< "Dumping Scotch graph file to " << str.name() << endl
                << "Use this in combination with dgpart." << endl;

            globalIndex globalCells(xadjSize-1);

            // Distributed graph file (.grf)
            label version = 2;
            str << version << nl;
            // Number of files (procglbnbr)
            str << Pstream::nProcs();
            // My file number (procloc)
            str << ' ' << Pstream::myProcNo() << nl;

            // Total number of vertices (vertglbnbr)
            str << globalCells.size();
            // Total number of connections (edgeglbnbr)
            str << ' ' << returnReduce(xadj[xadjSize-1], sumOp<label>())
                << nl;
            // Local number of vertices (vertlocnbr)
            str << xadjSize-1;
            // Local number of connections (edgelocnbr)
            str << ' ' << xadj[xadjSize-1] << nl;
            // Numbering starts from 0
            label baseval = 0;
            // 100*hasVertlabels+10*hasEdgeWeights+1*hasVertWeighs
            str << baseval << ' ' << "000" << nl;
            for (label celli = 0; celli < xadjSize-1; celli++)
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

    // Strategy
    // ~~~~~~~~

    // Default.
    SCOTCH_Strat stradat;
    check(SCOTCH_stratInit(&stradat), "SCOTCH_stratInit");

    if (decompositionDict_.found("scotchCoeffs"))
    {
        const dictionary& scotchCoeffs =
            decompositionDict_.subDict("scotchCoeffs");

        string strategy;
        if (scotchCoeffs.readIfPresent("strategy", strategy))
        {
            if (debug)
            {
                Info<< "ptscotch : Using strategy " << strategy << endl;
            }
            SCOTCH_stratDgraphMap(&stradat, strategy.c_str());
        }
    }


    // Graph
    // ~~~~~

    labelList velotab;

    // Check for externally provided cellweights and if so initialise weights
    if (returnReduce(cellWeights.size(), sumOp<label>()))
    {
        if (cellWeights.size() != xadjSize-1)
        {
            FatalErrorInFunction
                << "Number of cell weights " << cellWeights.size()
                << " does not equal number of cells " << xadjSize-1
                << exit(FatalError);
        }

        label nWeights = 1;
        velotab = scaleWeights(cellWeights, nWeights);

        if (nWeights == 1 && !cellWeights.size())
        {
            // Locally zero cells but not globally. Make sure we have
            // some size so .begin() does not return null pointer. Data
            // itself is never used.
            velotab.setSize(1, 1);
        }
    }


    if (debug)
    {
        Pout<< "SCOTCH_dgraphInit" << endl;
    }

    SCOTCH_Dgraph grafdat;
    check
    (
        SCOTCH_dgraphInit(&grafdat, PstreamGlobals::MPI_COMM_FOAM),
        "SCOTCH_dgraphInit"
    );


    if (debug)
    {
        Pout<< "SCOTCH_dgraphBuild with:" << nl
            << "xadjSize-1      : " << xadjSize-1 << nl
            << "xadj            : " << uintptr_t(xadj) << nl
            << "velotab         : " << uintptr_t(velotab.begin()) << nl
            << "adjncySize      : " << adjncySize << nl
            << "adjncy          : " << uintptr_t(adjncy) << nl
            << endl;
    }

    check
    (
        SCOTCH_dgraphBuild
        (
            &grafdat,               // grafdat
            0,                      // baseval, c-style numbering
            xadjSize-1,             // vertlocnbr, nCells
            xadjSize-1,             // vertlocmax
            const_cast<SCOTCH_Num*>(xadj),
                                    // vertloctab, start index per cell into
                                    // adjncy
            const_cast<SCOTCH_Num*>(xadj+1),// vendloctab, end index  ,,

            const_cast<SCOTCH_Num*>(velotab.begin()),// veloloctab, vtx weights
            nullptr,                   // vlblloctab

            adjncySize,             // edgelocnbr, number of arcs
            adjncySize,             // edgelocsiz
            const_cast<SCOTCH_Num*>(adjncy),         // edgeloctab
            nullptr,                   // edgegsttab
            nullptr                    // edlotab, edge weights
        ),
        "SCOTCH_dgraphBuild"
    );


    if (debug)
    {
        Pout<< "SCOTCH_dgraphCheck" << endl;
    }
    check(SCOTCH_dgraphCheck(&grafdat), "SCOTCH_dgraphCheck");


    // Architecture
    // ~~~~~~~~~~~~
    // (fully connected network topology since using switch)

    if (debug)
    {
        Pout<< "SCOTCH_archInit" << endl;
    }

    SCOTCH_Arch archdat;
    check(SCOTCH_archInit(&archdat), "SCOTCH_archInit");

    labelList processorWeights;
    if (decompositionDict_.found("scotchCoeffs"))
    {
        const dictionary& scotchCoeffs =
            decompositionDict_.subDict("scotchCoeffs");

        scotchCoeffs.readIfPresent("processorWeights", processorWeights);
    }
    if (processorWeights.size())
    {
        if (debug)
        {
            Info<< "ptscotch : Using processor weights "
                << processorWeights
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
        if (debug)
        {
            Pout<< "SCOTCH_archCmplt" << endl;
        }
        check
        (
            SCOTCH_archCmplt(&archdat, nProcessors_),
            "SCOTCH_archCmplt"
        );
    }


    // Hack:switch off fpu error trapping
    #ifdef  FE_NOMASK_ENV
    int oldExcepts = fedisableexcept
    (
        FE_DIVBYZERO
      | FE_INVALID
      | FE_OVERFLOW
    );
    #endif


    // Note: always provide allocated storage even if local size 0
    decomp.setSize(max(1, xadjSize-1));
    decomp = 0;

    if (debug)
    {
        Pout<< "SCOTCH_dgraphMap" << endl;
    }
    check
    (
        SCOTCH_dgraphMap
        (
            &grafdat,
            &archdat,
            &stradat,           // const SCOTCH_Strat *
            decomp.begin() // parttab
        ),
        "SCOTCH_graphMap"
    );

    #ifdef  FE_NOMASK_ENV
    feenableexcept(oldExcepts);
    #endif

    if (debug)
    {
        Pout<< "SCOTCH_dgraphExit" << endl;
    }

    // Release storage for graph
    SCOTCH_dgraphExit(&grafdat);

    // Release storage for strategy
    SCOTCH_stratExit(&stradat);

    // Release storage for network topology
    SCOTCH_archExit(&archdat);

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionMethods::ptscotch::ptscotch
(
    const dictionary& decompositionDict
)
:
    decompositionMethod(decompositionDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::decompositionMethods::ptscotch::decompose
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


Foam::labelList Foam::decompositionMethods::ptscotch::decompose
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

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
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

    // Rework back into decomposition for original mesh
    labelList fineDistribution(agglom.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = decomp[agglom[i]];
    }

    return fineDistribution;
}


Foam::labelList Foam::decompositionMethods::ptscotch::decompose
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
        "ptscotch",
        cellCells.m(),
        cellCells.offsets(),
        cellWeights,
        decomp
    );

    return decomp;
}


// ************************************************************************* //
