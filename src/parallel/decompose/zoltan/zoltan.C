/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2025 OpenFOAM Foundation
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

#include "zoltan.H"
#include "globalIndex.H"
#include "PstreamGlobals.H"
#include "Tuple3.H"
#include "addToRunTimeSelectionTable.H"

#include "sigFpe.H"
#ifdef LINUX_GNUC
    #ifndef __USE_GNU
        #define __USE_GNU
    #endif
    #include <fenv.h>
#endif

#include "zoltan.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionMethods
{
    defineTypeNameAndDebug(zoltan, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        zoltan,
        distributor
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

static int get_number_of_vertices(void *data, int *ierr)
{
    const Foam::pointField& points =
        *static_cast<const Foam::pointField*>(data);

    *ierr = ZOLTAN_OK;

    return points.size();
}


static void get_vertex_list
(
    void* data,
    int nGID,
    int nLID,
    ZOLTAN_ID_PTR globalIDs,
    ZOLTAN_ID_PTR localIDs,
    int wgt_dim,
    float* obj_wgts,
    int* ierr
)
{
    const Foam::Tuple3
    <
        const Foam::pointField&,
        const Foam::scalarField&,
        const Foam::globalIndex&
    >& vertexData =
        *static_cast
        <
            const Foam::Tuple3
            <
                const Foam::pointField&,
                const Foam::scalarField&,
                const Foam::globalIndex&
            >*
        >(data);

    const Foam::pointField& points = vertexData.first();
    const Foam::scalarField& weights = vertexData.second();
    const Foam::globalIndex& globalMap = vertexData.third();

    if (points.size() && wgt_dim != weights.size()/points.size())
    {
        *ierr = ZOLTAN_FATAL;
        return;
    }

    for (Foam::label i=0; i<points.size(); i++)
    {
        localIDs[i] = i;
        globalIDs[i] = globalMap.toGlobal(i);

        for(int j=0; j<wgt_dim; j++)
        {
            obj_wgts[wgt_dim*i + j] = weights[wgt_dim*i + j];
        }
    }

    *ierr = ZOLTAN_OK;
}


static int get_mesh_dim(void* data, int* ierr)
{
    return 3;
}


static void get_geom_list
(
    void* data,
    int nGID,
    int nLID,
    int nPoints,
    ZOLTAN_ID_PTR globalIDs,
    ZOLTAN_ID_PTR localIDs,
    int nDim,
    double* vertices,
    int* ierr
)
{
    const Foam::pointField& points =
        *static_cast<const Foam::pointField*>(data);

    if
    (
        (nGID != 1)
     || (nLID != 1)
     || (nPoints != points.size())
     || (nDim != 3)
    )
    {
        *ierr = ZOLTAN_FATAL;
        return;
    }

    double* p = vertices;

    for (Foam::label celli=0; celli<nPoints; celli++)
    {
        const Foam::point& pt = points[celli];

        for (Foam::direction cmpt=0; cmpt<Foam::vector::nComponents; cmpt++)
        {
            *p++ = pt[cmpt];
        }
    }

    *ierr = ZOLTAN_OK;
}


static void get_num_edges_list
(
    void* data,
    int nGID,
    int nLID,
    int nPoints,
    ZOLTAN_ID_PTR globalIDs,
    ZOLTAN_ID_PTR localIDs,
    int* numEdges,
    int* ierr
)
{
    const Foam::CompactListList<Foam::label>& adjacency =
        *static_cast<const Foam::CompactListList<Foam::label>*>(data);

    Foam::labelList sizes(adjacency.sizes());

    if ((nGID != 1) || (nLID != 1) || (nPoints != sizes.size()))
    {
        *ierr = ZOLTAN_FATAL;
        return;
    }

    for (Foam::label i=0; i<nPoints; i++)
    {
        numEdges[i] = sizes[i];
    }

    *ierr = ZOLTAN_OK;
}


static void get_edge_list
(
    void* data,
    int nGID,
    int nLID,
    int nPoints,
    ZOLTAN_ID_PTR globalIDs,
    ZOLTAN_ID_PTR localIDs,
    int* num_edges,
    ZOLTAN_ID_PTR nborGID,
    int* nborProc,
    int wgt_dim,
    float* ewgts,
    int* ierr
)
{
    const Foam::Tuple2
    <
        const Foam::CompactListList<Foam::label>&,
        const Foam::globalIndex&
    >& edgeData =
        *static_cast
        <
            const Foam::Tuple2
            <
                const Foam::CompactListList<Foam::label>&,
                const Foam::globalIndex&
            >*
        >(data);

    const Foam::CompactListList<Foam::label>& adjacency = edgeData.first();
    const Foam::globalIndex& globalMap = edgeData.second();

    const Foam::labelList& m(adjacency.m());

    if
    (
        (nGID != 1)
     || (nLID != 1)
    )
    {
        *ierr = ZOLTAN_FATAL;
        return;
    }

    forAll(m, i)
    {
        nborGID[i] = m[i];
        nborProc[i] = globalMap.whichProcID(m[i]);
    }

    *ierr = ZOLTAN_OK;
}


Foam::label Foam::decompositionMethods::zoltan::decompose
(
    const CompactListList<label>& adjacency,
    const pointField& points,
    const scalarField& pWeights,
    List<label>& decomp
) const
{
    stringList args(1);
    args[0] = "zoltan";

    const int argc = args.size();
    List<char*> argv(argc);
    for (label i = 0; i < argc; i++)
    {
        argv[i] = strdup(args[i].c_str());
    }

    float ver;
    int rc = Zoltan_Initialize(argc, argv.begin(), &ver);

    if (rc != ZOLTAN_OK)
    {
        FatalErrorInFunction
            << "Failed initialising Zoltan" << exit(FatalError);
    }

    struct Zoltan_Struct *zz = Zoltan_Create(PstreamGlobals::MPI_COMM_FOAM);

    const int nWeights = this->nWeights(points, pWeights);

    // Set internal parameters
    Zoltan_Set_Param(zz, "return_lists", "export");
    Zoltan_Set_Param(zz, "obj_weight_dim", name(nWeights).c_str());
    Zoltan_Set_Param(zz, "edge_weight_dim", "0");

    // General default parameters
    Zoltan_Set_Param(zz, "debug_level", "0");
    Zoltan_Set_Param(zz, "imbalance_tol", "1.05");

    word lb_method("graph");

    // Set default method parameters
    Zoltan_Set_Param(zz, "lb_method", lb_method.c_str());
    Zoltan_Set_Param(zz, "lb_approach", "repartition");

    if (!methodDict_.empty())
    {
        methodDict_.readIfPresent("lb_method", lb_method);

        forAllConstIter(IDLList<entry>, methodDict_, iter)
        {
            if (!iter().isDict())
            {
                const word& key = iter().keyword();
                const word value(iter().stream());

                if (debug)
                {
                    Info<< typeName << " : setting parameter " << key
                        << " to " << value << endl;
                }

                Zoltan_Set_Param(zz, key.c_str(), value.c_str());
            }
        }
    }

    if (nWeights > 1 && lb_method != "rcb")
    {
        FatalIOErrorInFunction(methodDict_)
            << "Multiple constraints specified for lb_method " << lb_method
            << " but is only supported by the rcb method"
            << exit(FatalIOError);
    }

    globalIndex globalMap(points.size());

    void* pointsPtr = &const_cast<pointField&>(points);
    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, pointsPtr);

    Tuple3<const pointField&, const scalarField&, const globalIndex&> vertexData
    (
        points,
        pWeights,
        globalMap
    );
    Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, &vertexData);

    // Callbacks for geometry
    Zoltan_Set_Num_Geom_Fn(zz, get_mesh_dim, pointsPtr);
    Zoltan_Set_Geom_Multi_Fn(zz, get_geom_list, pointsPtr);

    // Callbacks for connectivity
    void* adjacencyPtr = &const_cast<CompactListList<label>&>(adjacency);
    Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list, adjacencyPtr);

    Tuple2<const CompactListList<label>&, const globalIndex&> edgeData
    (
        adjacency,
        globalMap
    );
    Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list, &edgeData);

    int changed, num_imp, num_exp, *imp_procs, *exp_procs;
    int *imp_to_part, *exp_to_part;
    int num_gid_entries, num_lid_entries;
    ZOLTAN_ID_PTR imp_global_ids, exp_global_ids;
    ZOLTAN_ID_PTR imp_local_ids, exp_local_ids;

    // Switch off FPU error trapping to work around divide-by-zero bug
    // in the Zoltan graph and hypergraph methods
    #ifdef FE_NOMASK_ENV
    int oldExcepts = fedisableexcept
    (
        FE_DIVBYZERO
      | FE_INVALID
      | FE_OVERFLOW
    );
    #endif

    // Perform load balancing
    Zoltan_LB_Partition
    (
        zz,
        &changed,
        &num_gid_entries,
        &num_lid_entries,
        &num_imp,
        &imp_global_ids,
        &imp_local_ids,
        &imp_procs,
        &imp_to_part,
        &num_exp,
        &exp_global_ids,
        &exp_local_ids,
        &exp_procs,
        &exp_to_part
    );

    // Reinstate FPU error trapping
    #ifdef FE_NOMASK_ENV
    feenableexcept(oldExcepts);
    #endif

    if (debug && changed)
    {
        Pout << "Zoltan number to move " << changed << " " << num_exp << endl;
    }

    decomp.setSize(points.size());
    decomp = Pstream::myProcNo();

    if (changed)
    {
        for(int i=0; i<num_exp; i++)
        {
            decomp[exp_local_ids[i]] = exp_procs[i];
        }

        // Sum the number of cells allocated to each processor
        labelList nProcCells(Pstream::nProcs(), 0);

        forAll(decomp, i)
        {
            nProcCells[decomp[i]]++;
        }

        reduce(nProcCells, ListOp<sumOp<label>>());

        // If there are no cells allocated to this processor keep the first one
        // to ensure that all processors have at least one cell
        if (nProcCells[Pstream::myProcNo()] == 0)
        {
            Pout<< "    No cells allocated to this processor"
                   ", keeping first cell"
                << endl;
            decomp[0] = Pstream::myProcNo();
        }
    }

    // Free memory allocated for load-balancing results by Zoltan library

    Zoltan_LB_Free_Part
    (
        &imp_global_ids,
        &imp_local_ids,
        &imp_procs,
        &imp_to_part
    );

    Zoltan_LB_Free_Part
    (
        &exp_global_ids,
        &exp_local_ids,
        &exp_procs,
        &exp_to_part
    );

    Zoltan_Destroy(&zz);

    // Free the string storage allocated by strdup
    for (label i = 0; i < argc; i++)
    {
        free(argv[i]);
    }

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionMethods::zoltan::zoltan
(
    const dictionary& decompositionDict,
    const dictionary& methodDict
)
:
    decompositionMethod(decompositionDict),
    methodDict_(methodDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::decompositionMethods::zoltan::decompose
(
    const polyMesh& mesh,
    const pointField& cellCentres,
    const scalarField& cellWeights
)
{
    if (cellCentres.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Can use this decomposition method only for the whole mesh"
            << endl
            << "and supply one coordinate (cellCentre) for every cell." << endl
            << "The number of coordinates " << cellCentres.size() << endl
            << "The number of cells in the mesh " << mesh.nCells()
            << exit(FatalError);
    }


    // Convert mesh.cellCells() into compressed format
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
    List<label> decomp;
    decompose
    (
        cellCells,
        cellCentres,
        cellWeights,
        decomp
    );

    return decomp;
}


Foam::labelList Foam::decompositionMethods::zoltan::decompose
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

    // Convert agglom into compressed format
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
    List<label> decomp;
    decompose
    (
        cellCells,
        agglomPoints,
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


Foam::labelList Foam::decompositionMethods::zoltan::decompose
(
    const labelListList& globalPointPoints,
    const pointField& points,
    const scalarField& pWeights
)
{
    if (points.size() != globalPointPoints.size())
    {
        FatalErrorInFunction
            << "Inconsistent number of cells (" << globalPointPoints.size()
            << ") and number of cell centres (" << points.size()
            << ")." << exit(FatalError);
    }


    // Convert globalPointPoints into compressed format
    CompactListList<label> pointPoints(globalPointPoints);

    // Decompose using weights
    List<label> decomp;
    decompose
    (
        pointPoints,
        points,
        pWeights,
        decomp
    );

    return decomp;
}


// ************************************************************************* //
