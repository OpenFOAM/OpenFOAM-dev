/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "metisDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "floatScalar.H"
#include "Time.H"

extern "C"
{
#define OMPI_SKIP_MPICXX
#   include "metis.h"
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(metisDecomp, 0);

    addToRunTimeSelectionTable(decompositionMethod, metisDecomp, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Call Metis with options from dictionary.
Foam::label Foam::metisDecomp::decompose
(
    const List<int>& adjncy,
    const List<int>& xadj,
    const scalarField& cWeights,

    List<int>& finalDecomp
)
{
    // C style numbering
    //int numFlag = 0;

    // Method of decomposition
    // recursive: multi-level recursive bisection (default)
    // k-way: multi-level k-way
    word method("recursive");

    int numCells = xadj.size()-1;

    // decomposition options
    List<int> options(METIS_NOPTIONS);
    METIS_SetDefaultOptions(options.begin());

    // processor weights initialised with no size, only used if specified in
    // a file
    Field<floatScalar> processorWeights;

    // cell weights (so on the vertices of the dual)
    List<int> cellWeights;

    // face weights (so on the edges of the dual)
    List<int> faceWeights;


    // Check for externally provided cellweights and if so initialise weights
    scalar minWeights = gMin(cWeights);
    if (cWeights.size() > 0)
    {
        if (minWeights <= 0)
        {
            WarningIn
            (
                "metisDecomp::decompose"
                "(const pointField&, const scalarField&)"
            )   << "Illegal minimum weight " << minWeights
                << endl;
        }

        if (cWeights.size() != numCells)
        {
            FatalErrorIn
            (
                "metisDecomp::decompose"
                "(const pointField&, const scalarField&)"
            )   << "Number of cell weights " << cWeights.size()
                << " does not equal number of cells " << numCells
                << exit(FatalError);
        }
        // Convert to integers.
        cellWeights.setSize(cWeights.size());
        forAll(cellWeights, i)
        {
            cellWeights[i] = int(cWeights[i]/minWeights);
        }
    }


    // Check for user supplied weights and decomp options
    if (decompositionDict_.found("metisCoeffs"))
    {
        const dictionary& metisCoeffs =
            decompositionDict_.subDict("metisCoeffs");
        word weightsFile;

        if (metisCoeffs.readIfPresent("method", method))
        {
            if (method != "recursive" && method != "k-way")
            {
                FatalErrorIn("metisDecomp::decompose()")
                    << "Method " << method << " in metisCoeffs in dictionary : "
                    << decompositionDict_.name()
                    << " should be 'recursive' or 'k-way'"
                    << exit(FatalError);
            }

            Info<< "metisDecomp : Using Metis method     " << method
                << nl << endl;
        }

        if (metisCoeffs.readIfPresent("options", options))
        {
            if (options.size() != METIS_NOPTIONS)
            {
                FatalErrorIn("metisDecomp::decompose()")
                    << "Number of options in metisCoeffs in dictionary : "
                    << decompositionDict_.name()
                    << " should be " << METIS_NOPTIONS
                    << exit(FatalError);
            }

            Info<< "metisDecomp : Using Metis options     " << options
                << nl << endl;
        }

        if (metisCoeffs.readIfPresent("processorWeights", processorWeights))
        {
            processorWeights /= sum(processorWeights);

            if (processorWeights.size() != nProcessors_)
            {
                FatalErrorIn("metisDecomp::decompose(const pointField&)")
                    << "Number of processor weights "
                    << processorWeights.size()
                    << " does not equal number of domains " << nProcessors_
                    << exit(FatalError);
            }
        }

        //if (metisCoeffs.readIfPresent("cellWeightsFile", weightsFile))
        //{
        //    Info<< "metisDecomp : Using cell-based weights." << endl;
        //
        //    IOList<int> cellIOWeights
        //    (
        //        IOobject
        //        (
        //            weightsFile,
        //            mesh_.time().timeName(),
        //            mesh_,
        //            IOobject::MUST_READ,
        //            IOobject::AUTO_WRITE
        //        )
        //    );
        //    cellWeights.transfer(cellIOWeights);
        //
        //    if (cellWeights.size() != xadj.size()-1)
        //    {
        //        FatalErrorIn("metisDecomp::decompose(const pointField&)")
        //            << "Number of cell weights " << cellWeights.size()
        //            << " does not equal number of cells " << xadj.size()-1
        //            << exit(FatalError);
        //    }
        //}
    }

    int ncon = 1;

    int nProcs = nProcessors_;

    // output: cell -> processor addressing
    finalDecomp.setSize(numCells);

    // output: number of cut edges
    int edgeCut = 0;

    if (method == "recursive")
    {
        METIS_PartGraphRecursive
        (
            &numCells,          // num vertices in graph
            &ncon,              // num balancing constraints
            const_cast<List<int>&>(xadj).begin(),   // indexing into adjncy
            const_cast<List<int>&>(adjncy).begin(), // neighbour info
            cellWeights.begin(),// vertexweights
            NULL,               // vsize: total communication vol
            faceWeights.begin(),// edgeweights
            &nProcs,            // nParts
            processorWeights.begin(),   // tpwgts
            NULL,               // ubvec: processor imbalance (default)
            options.begin(),
            &edgeCut,
            finalDecomp.begin()
        );
    }
    else
    {
        METIS_PartGraphKway
        (
            &numCells,         // num vertices in graph
            &ncon,              // num balancing constraints
            const_cast<List<int>&>(xadj).begin(),   // indexing into adjncy
            const_cast<List<int>&>(adjncy).begin(), // neighbour info
            cellWeights.begin(),// vertexweights
            NULL,               // vsize: total communication vol
            faceWeights.begin(),// edgeweights
            &nProcs,            // nParts
            processorWeights.begin(),   // tpwgts
            NULL,               // ubvec: processor imbalance (default)
            options.begin(),
            &edgeCut,
            finalDecomp.begin()
        );
    }

    return edgeCut;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::metisDecomp::metisDecomp(const dictionary& decompositionDict)
:
    decompositionMethod(decompositionDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::metisDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
)
{
    if (points.size() != mesh.nCells())
    {
        FatalErrorIn
        (
            "metisDecomp::decompose(const pointField&,const scalarField&)"
        )   << "Can use this decomposition method only for the whole mesh"
            << endl
            << "and supply one coordinate (cellCentre) for every cell." << endl
            << "The number of coordinates " << points.size() << endl
            << "The number of cells in the mesh " << mesh.nCells()
            << exit(FatalError);
    }

    CompactListList<label> cellCells;
    calcCellCells
    (
        mesh,
        identity(mesh.nCells()),
        mesh.nCells(),
        false,
        cellCells
    );

    // Decompose using default weights
    labelList decomp;
    decompose(cellCells.m(), cellCells.offsets(), pointWeights, decomp);

    return decomp;
}


Foam::labelList Foam::metisDecomp::decompose
(
    const polyMesh& mesh,
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& agglomWeights
)
{
    if (agglom.size() != mesh.nCells())
    {
        FatalErrorIn
        (
            "metisDecomp::decompose"
            "(const labelList&, const pointField&, const scalarField&)"
        )   << "Size of cell-to-coarse map " << agglom.size()
            << " differs from number of cells in mesh " << mesh.nCells()
            << exit(FatalError);
    }

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    CompactListList<label> cellCells;
    calcCellCells(mesh, agglom, agglomPoints.size(), false, cellCells);

    // Decompose using default weights
    labelList finalDecomp;
    decompose(cellCells.m(), cellCells.offsets(), agglomWeights, finalDecomp);


    // Rework back into decomposition for original mesh
    labelList fineDistribution(agglom.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = finalDecomp[agglom[i]];
    }

    return finalDecomp;
}


Foam::labelList Foam::metisDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cellWeights
)
{
    if (cellCentres.size() != globalCellCells.size())
    {
        FatalErrorIn
        (
            "metisDecomp::decompose"
            "(const pointField&, const labelListList&, const scalarField&)"
        )   << "Inconsistent number of cells (" << globalCellCells.size()
            << ") and number of cell centres (" << cellCentres.size()
            << ")." << exit(FatalError);
    }


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
