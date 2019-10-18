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

#include "scotchDecomp.H"
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
    defineTypeNameAndDebug(scotchDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        scotchDecomp,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::scotchDecomp::check(const int retVal, const char* str)
{
    if (retVal)
    {
        FatalErrorInFunction
            << "Call to scotch routine " << str << " failed."
            << exit(FatalError);
    }
}


Foam::label Foam::scotchDecomp::decompose
(
    const fileName& meshPath,
    const List<label>& adjncy,
    const List<label>& xadj,
    const scalarField& cWeights,

    List<label>& finalDecomp
)
{
    if (!Pstream::parRun())
    {
        decomposeOneProc
        (
            meshPath,
            adjncy,
            xadj,
            cWeights,
            finalDecomp
        );
    }
    else
    {
        if (debug)
        {
            Info<< "scotchDecomp : running in parallel."
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
            forAll(cWeights, celli)
            {
                allXadj[nTotalCells] = xadj[celli];
                allWeights[nTotalCells++] = cWeights[celli];
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
            finalDecomp = SubField<label>
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
                    << cWeights;
            }

            // Receive back decomposition
            IPstream fromMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );
            fromMaster >> finalDecomp;
        }
    }
    return 0;
}


// Call scotch with options from dictionary.
Foam::label Foam::scotchDecomp::decomposeOneProc
(
    const fileName& meshPath,
    const List<label>& adjncy,
    const List<label>& xadj,
    const scalarField& cWeights,

    List<label>& finalDecomp
)
{
    // Dump graph
    if (decompositionDict_.found("scotchCoeffs"))
    {
        const dictionary& scotchCoeffs =
            decompositionDict_.subDict("scotchCoeffs");

        if (scotchCoeffs.lookupOrDefault("writeGraph", false))
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
                Info<< "scotchDecomp : Using strategy " << strategy << endl;
            }
            SCOTCH_stratGraphMap(&stradat, strategy.c_str());
            // fprintf(stdout, "S\tStrat=");
            // SCOTCH_stratSave(&stradat, stdout);
            // fprintf(stdout, "\n");
        }
    }


    // Graph
    // ~~~~~

    List<label> velotab;


    // Check for externally provided cellweights and if so initialise weights
    // Note: min, not gMin since routine runs on master only.
    scalar minWeights = min(cWeights);
    if (!cWeights.empty())
    {
        if (minWeights <= 0)
        {
            WarningInFunction
                << "Illegal minimum weight " << minWeights
                << endl;
        }

        if (cWeights.size() != xadj.size()-1)
        {
            FatalErrorInFunction
                << "Number of cell weights " << cWeights.size()
                << " does not equal number of cells " << xadj.size()-1
                << exit(FatalError);
        }

        scalar velotabSum = sum(cWeights)/minWeights;

        scalar rangeScale(1.0);

        if (velotabSum > scalar(labelMax - 1))
        {
            // 0.9 factor of safety to avoid floating point round-off in
            // rangeScale tipping the subsequent sum over the integer limit.
            rangeScale = 0.9*scalar(labelMax - 1)/velotabSum;

            WarningInFunction
                << "Sum of weights has overflowed integer: " << velotabSum
                << ", compressing weight scale by a factor of " << rangeScale
                << endl;
        }

        // Convert to integers.
        velotab.setSize(cWeights.size());

        forAll(velotab, i)
        {
            velotab[i] = int((cWeights[i]/minWeights - 1)*rangeScale) + 1;
        }
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
            nullptr,                   // vlbltab
            adjncy.size(),          // edgenbr, number of arcs
            adjncy.begin(),         // edgetab
            nullptr                    // edlotab, edge weights
        ),
        "SCOTCH_graphBuild"
    );
    check(SCOTCH_graphCheck(&grafdat), "SCOTCH_graphCheck");


    // Architecture
    // ~~~~~~~~~~~~
    // (fully connected network topology since using switch)

    SCOTCH_Arch archdat;
    check(SCOTCH_archInit(&archdat), "SCOTCH_archInit");

    List<label> processorWeights;
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
            Info<< "scotchDecomp : Using processor weights " << processorWeights
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


        //- Hack to test clustering. Note that finalDecomp is non-compact
        //  numbers!
        //
        ////- Set up variable sizes architecture
        // check
        //(
        //    SCOTCH_archVcmplt(&archdat),
        //    "SCOTCH_archVcmplt"
        //);
        //
        ////- Stategy flags: go for quality or load balance (or leave default)
        // SCOTCH_Num straval = 0;
        ////straval |= SCOTCH_STRATQUALITY;
        ////straval |= SCOTCH_STRATQUALITY;
        //
        ////- Number of cells per agglomeration
        ////SCOTCH_Num agglomSize = SCOTCH_archSize(&archdat);
        // SCOTCH_Num agglomSize = 3;
        //
        ////- Build strategy for agglomeration
        // check
        //(
        //    SCOTCH_stratGraphClusterBuild
        //    (
        //        &stradat,   // strategy to build
        //        straval,    // strategy flags
        //        agglomSize, // cells per cluster
        //        1.0,        // weight?
        //        0.01        // max load imbalance
        //    ),
        //    "SCOTCH_stratGraphClusterBuild"
        //);
    }


    // SCOTCH_Mapping mapdat;
    // SCOTCH_graphMapInit(&grafdat, &mapdat, &archdat, nullptr);
    // SCOTCH_graphMapCompute(&grafdat, &mapdat, &stradat); // Perform mapping
    // SCOTCH_graphMapExit(&grafdat, &mapdat);


    // Hack:switch off fpu error trapping
    #ifdef FE_NOMASK_ENV
    int oldExcepts = fedisableexcept
    (
        FE_DIVBYZERO
      | FE_INVALID
      | FE_OVERFLOW
    );
    #endif

    finalDecomp.setSize(xadj.size()-1);
    finalDecomp = 0;
    check
    (
        SCOTCH_graphMap
        (
            &grafdat,
            &archdat,
            &stradat,           // const SCOTCH_Strat *
            finalDecomp.begin() // parttab
        ),
        "SCOTCH_graphMap"
    );

    #ifdef FE_NOMASK_ENV
    feenableexcept(oldExcepts);
    #endif



    // finalDecomp.setSize(xadj.size()-1);
    // check
    //(
    //    SCOTCH_graphPart
    //    (
    //        &grafdat,
    //        nProcessors_,       // partnbr
    //        &stradat,           // const SCOTCH_Strat *
    //        finalDecomp.begin() // parttab
    //    ),
    //    "SCOTCH_graphPart"
    //);

    // Release storage for graph
    SCOTCH_graphExit(&grafdat);
    // Release storage for strategy
    SCOTCH_stratExit(&stradat);
    // Release storage for network topology
    SCOTCH_archExit(&archdat);

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scotchDecomp::scotchDecomp(const dictionary& decompositionDict)
:
    decompositionMethod(decompositionDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::scotchDecomp::decompose
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

    // Calculate local or global (if Pstream::parRun()) connectivity
    CompactListList<label> cellCells;
    calcCellCells
    (
        mesh,
        identity(mesh.nCells()),
        mesh.nCells(),
        true,
        cellCells
    );

    // Decompose using default weights
    List<label> finalDecomp;
    decompose
    (
        mesh.time().path()/mesh.name(),
        cellCells.m(),
        cellCells.offsets(),
        pointWeights,
        finalDecomp
    );

    // Copy back to labelList
    labelList decomp(finalDecomp.size());
    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }
    return decomp;
}


Foam::labelList Foam::scotchDecomp::decompose
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
    List<label> finalDecomp;
    decompose
    (
        mesh.time().path()/mesh.name(),
        cellCells.m(),
        cellCells.offsets(),
        pointWeights,
        finalDecomp
    );

    // Rework back into decomposition for original mesh_
    labelList fineDistribution(agglom.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = finalDecomp[agglom[i]];
    }

    return fineDistribution;
}


Foam::labelList Foam::scotchDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cWeights
)
{
    if (cellCentres.size() != globalCellCells.size())
    {
        FatalErrorInFunction
            << "Inconsistent number of cells (" << globalCellCells.size()
            << ") and number of cell centres (" << cellCentres.size()
            << ")." << exit(FatalError);
    }


    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    CompactListList<label> cellCells(globalCellCells);

    // Decompose using weights
    List<label> finalDecomp;
    decompose
    (
        "scotch",
        cellCells.m(),
        cellCells.offsets(),
        cWeights,
        finalDecomp
    );

    // Copy back to labelList
    labelList decomp(finalDecomp.size());
    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }
    return decomp;
}


// ************************************************************************* //
