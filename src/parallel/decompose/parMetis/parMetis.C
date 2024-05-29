/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "parMetis.H"
#include "Time.H"
#include "globalIndex.H"
#include "labelIOField.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionMethods
{
    defineTypeNameAndDebug(parMetis, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        parMetis,
        distributor
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::decompositionMethods::parMetis::decompose
(
    const labelList& xadj,
    const labelList& adjncy,
    const pointField& cellCentres,
    const label nWeights,
    const labelList& cellWeights,
    const labelList& faceWeights,
    labelList& decomp
)
{
    // C style numbering
    label numFlag = 0;

    List<real_t> processorWeights;

    if (processorWeights_.size())
    {
        if (processorWeights_.size() != nWeights*nProcessors_)
        {
            FatalErrorInFunction
                << "Number of processor weights specified in parMetisCoeffs "
                << processorWeights_.size()
                << " does not equal number of constraints * number of domains "
                << nWeights*nProcessors_
                << exit(FatalError);
        }

        processorWeights = processorWeights_;
    }
    else
    {
        processorWeights.setSize(nWeights*nProcessors_, 1.0/nProcessors_);
    }

    // Imbalance tolerance
    List<real_t> ubvec(nWeights, 1.02);

    // If only one processor there is no imbalance
    if (nProcessors_ == 1)
    {
        ubvec[0] = 1;
    }

    // Distribute to all processors the number of cells on each processor
    globalIndex globalMap(cellCentres.size());

    // Get the processor-cell offset table
    labelList cellOffsets(globalMap.offsets());

    // Weight info
    label wgtFlag = 0;
    const label* vwgtPtr = nullptr;
    const label* adjwgtPtr = nullptr;

    // Weights on vertices of graph (cells)
    if (cellWeights.size())
    {
        vwgtPtr = cellWeights.begin();
        wgtFlag += 2;
    }

    // Weights on edges of graph (faces)
    if (faceWeights.size())
    {
        adjwgtPtr = faceWeights.begin();
        wgtFlag += 1;
    }

    MPI_Comm comm = MPI_COMM_WORLD;

    // Output: cell -> processor addressing
    decomp.setSize(cellCentres.size());

    // Output: the number of edges that are cut by the partitioning
    label edgeCut = 0;

    if (method_ == "kWay")
    {
        ParMETIS_V3_PartKway
        (
            cellOffsets.begin(),
            const_cast<label*>(xadj.begin()),
            const_cast<label*>(adjncy.begin()),
            const_cast<label*>(vwgtPtr),
            const_cast<label*>(adjwgtPtr),
            &wgtFlag,
            &numFlag,
            const_cast<label*>(&nWeights),
            &nProcessors_,
            processorWeights.begin(),
            ubvec.begin(),
            const_cast<labelList&>(options_).begin(),
            &edgeCut,
            decomp.begin(),
            &comm
        );
    }
    else if (method_ == "geomKway")
    {
        // Number of dimensions
        label nDims = 3;

        // Convert pointField into List<real_t>
        List<real_t> xyz(nDims*cellCentres.size());
        label i = 0;
        forAll(cellCentres, celli)
        {
            const point& cc = cellCentres[celli];
            xyz[i++] = cc.x();
            xyz[i++] = cc.y();
            xyz[i++] = cc.z();
        }

        ParMETIS_V3_PartGeomKway
        (
            cellOffsets.begin(),
            const_cast<label*>(xadj.begin()),
            const_cast<label*>(adjncy.begin()),
            const_cast<label*>(vwgtPtr),
            const_cast<label*>(adjwgtPtr),
            &wgtFlag,
            &numFlag,
            &nDims,
            xyz.begin(),
            const_cast<label*>(&nWeights),
            &nProcessors_,
            processorWeights.begin(),
            ubvec.begin(),
            const_cast<labelList&>(options_).begin(),
            &edgeCut,
            decomp.begin(),
            &comm
        );
    }
    else if (method_ == "adaptiveRepart")
    {
        // Size of the vertices with respect to redistribution cost
        labelList vsize(cellCentres.size(), 1);

        ParMETIS_V3_AdaptiveRepart
        (
            cellOffsets.begin(),
            const_cast<label*>(xadj.begin()),
            const_cast<label*>(adjncy.begin()),
            const_cast<label*>(vwgtPtr),
            const_cast<label*>(vsize.begin()),
            const_cast<label*>(adjwgtPtr),
            &wgtFlag,
            &numFlag,
            const_cast<label*>(&nWeights),
            &nProcessors_,
            processorWeights.begin(),
            ubvec.begin(),
            &itr_,
            const_cast<labelList&>(options_).begin(),
            &edgeCut,
            decomp.begin(),
            &comm
        );
    }

    return edgeCut;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionMethods::parMetis::parMetis
(
    const dictionary& decompositionDict
)
:
    decompositionMethod(decompositionDict),
    method_("geomKway"),
    options_(4, 0),
    itr_(1000)
{
    // Check for user supplied weights and decomp options
    if (decompositionDict.found("parMetisCoeffs"))
    {
        const dictionary& parMetisCoeffs =
            decompositionDict.subDict("parMetisCoeffs");

        Info<< type() << ": reading coefficients:" << endl;

        if (parMetisCoeffs.readIfPresent("method", method_))
        {
            if
            (
                method_ != "kWay"
             && method_ != "geomKWay"
             && method_ != "adaptiveRepart"
            )
            {
                FatalIOErrorInFunction(parMetisCoeffs)
                    << "Method " << method_
                    << " in parMetisCoeffs in dictionary : "
                    << decompositionDict.name()
                    << " should be kWay, geomKWay or adaptiveRepart"
                    << exit(FatalIOError);
            }

            Info<< "    method: " << method_ << endl;
        }

        if
        (
            method_ == "adaptiveRepart"
         && parMetisCoeffs.readIfPresent("itr", itr_)
        )
        {
            Info<< "    itr: " << itr_ << endl;
        }

        if (parMetisCoeffs.readIfPresent("options", options_))
        {
            if (options_.size() != 4)
            {
                FatalIOErrorInFunction(parMetisCoeffs)
                    << "Number of options in parMetisCoeffs dictionary : "
                    << decompositionDict.name()
                    << " should be 4, found " << options_
                    << exit(FatalIOError);
            }

            Info<< "    options: " << options_ << endl;
        }

        parMetisCoeffs.readIfPresent("processorWeights_", processorWeights_);

        Info << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::decompositionMethods::parMetis::decompose
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

    if (Pstream::parRun() && !points.size())
    {
        Pout<< "No points on processor " << Pstream::myProcNo() << nl
            << "    ParMETIS cannot redistribute without points"
               " on all processors." << endl;
        FatalErrorInFunction << exit(FatalError);
    }

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

    const label nWeights = this->nWeights(points, pointWeights);

    labelList decomp;
    decompose
    (
        cellCells.offsets(),
        cellCells.m(),
        points,
        nWeights,
        scaleWeights(pointWeights, nWeights),
        labelList(),
        decomp
    );

    return decomp;
}


Foam::labelList Foam::decompositionMethods::parMetis::decompose
(
    const polyMesh& mesh,
    const labelList& cellToRegion,
    const pointField& regionPoints,
    const scalarField& pointWeights
)
{
    if (cellToRegion.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Size of cell-to-coarse map " << cellToRegion.size()
            << " differs from number of cells in mesh " << mesh.nCells()
            << exit(FatalError);
    }

    if (Pstream::parRun() && !regionPoints.size())
    {
        Pout<< "No points on processor " << Pstream::myProcNo() << nl
            << "    ParMETIS cannot redistribute without points"
               " on all processors." << endl;
        FatalErrorInFunction << exit(FatalError);
    }

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
    CompactListList<label> cellCells;
    calcCellCells
    (
        mesh,
        cellToRegion,
        regionPoints.size(),
        true,
        cellCells
    );

    const label nWeights = this->nWeights(regionPoints, pointWeights);

    // Decompose using weights
    labelList decomp;
    decompose
    (
        cellCells.m(),
        cellCells.offsets(),
        regionPoints,
        nWeights,
        scaleWeights(pointWeights, nWeights),
        labelList(),
        decomp
    );

    // Rework back into decomposition for original mesh
    labelList fineDistribution(cellToRegion.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = decomp[cellToRegion[i]];
    }

    return fineDistribution;
}


Foam::labelList Foam::decompositionMethods::parMetis::decompose
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

    if (Pstream::parRun() && !cellCentres.size())
    {
        Pout<< "No points on processor " << Pstream::myProcNo() << nl
            << "    ParMETIS cannot redistribute without points"
               " on all processors." << endl;
        FatalErrorInFunction << exit(FatalError);
    }

    // Make Metis Distributed CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
    CompactListList<label> cellCells(globalCellCells);

    const label nWeights = this->nWeights(cellCentres, cellWeights);

    labelList decomp;
    decompose
    (
        cellCells.offsets(),
        cellCells.m(),
        cellCentres,
        nWeights,
        scaleWeights(cellWeights, nWeights),
        labelList(),
        decomp
    );

    return decomp;
}


// ************************************************************************* //
