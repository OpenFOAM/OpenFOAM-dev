/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    #include "metis.h"
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(metisDecomp, 0);
    addToRunTimeSelectionTable(decompositionMethod, metisDecomp, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::metisDecomp::decompose
(
    const List<label>& adjncy,
    const List<label>& xadj,
    const scalarField& cWeights,

    List<label>& finalDecomp
)
{
    // Method of decomposition
    // recursive: multi-level recursive bisection (default)
    // k-way: multi-level k-way
    word method("recursive");

    label numCells = xadj.size()-1;

    // decomposition options
    List<label> options(METIS_NOPTIONS);
    METIS_SetDefaultOptions(options.begin());

    // processor weights initialised with no size, only used if specified in
    // a file
    Field<floatScalar> processorWeights;

    // cell weights (so on the vertices of the dual)
    List<label> cellWeights;

    // face weights (so on the edges of the dual)
    List<label> faceWeights;


    // Check for externally provided cellweights and if so initialise weights
    scalar minWeights = gMin(cWeights);
    if (cWeights.size() > 0)
    {
        if (minWeights <= 0)
        {
            WarningInFunction
                << "Illegal minimum weight " << minWeights
                << endl;
        }

        if (cWeights.size() != numCells)
        {
            FatalErrorInFunction
                << "Number of cell weights " << cWeights.size()
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
                FatalErrorInFunction
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
                FatalErrorInFunction
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
                FatalErrorInFunction
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
        //    IOList<label> cellIOWeights
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
        //        FatalErrorInFunction
        //            << "Number of cell weights " << cellWeights.size()
        //            << " does not equal number of cells " << xadj.size()-1
        //            << exit(FatalError);
        //    }
        //}
    }

    label ncon = 1;

    label nProcs = nProcessors_;

    // output: cell -> processor addressing
    finalDecomp.setSize(numCells);

    // output: number of cut edges
    label edgeCut = 0;

    if (method == "recursive")
    {
        METIS_PartGraphRecursive
        (
            &numCells,          // num vertices in graph
            &ncon,              // num balancing constraints
            const_cast<List<label>&>(xadj).begin(),   // indexing into adjncy
            const_cast<List<label>&>(adjncy).begin(), // neighbour info
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
            const_cast<List<label>&>(xadj).begin(),   // indexing into adjncy
            const_cast<List<label>&>(adjncy).begin(), // neighbour info
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
        FatalErrorInFunction
            << "Can use this decomposition method only for the whole mesh"
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
        FatalErrorInFunction
            << "Size of cell-to-coarse map " << agglom.size()
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
        FatalErrorInFunction
            << "Inconsistent number of cells (" << globalCellCells.size()
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
