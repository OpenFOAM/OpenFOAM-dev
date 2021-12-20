/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "zoltanRenumber.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "labelIOList.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "globalIndex.H"
#include "uint.H"
#include "PstreamGlobals.H"

#include "zoltan.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoltanRenumber, 0);

    addToRunTimeSelectionTable
    (
        renumberMethod,
        zoltanRenumber,
        dictionary
    );
}


static int get_number_of_vertices(void *data, int *ierr)
{
    const Foam::polyMesh& mesh = *static_cast<const Foam::polyMesh*>(data);

    *ierr = ZOLTAN_OK;

    return mesh.nCells();
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
    const Foam::polyMesh& mesh = *static_cast<const Foam::polyMesh*>(data);

    // In this example, return the IDs of our vertices, but no weights.
    // Zoltan will assume equally weighted vertices.

    // Global calculation engine
    Foam::globalIndex globalCellMap(mesh.nCells());

    for (Foam::label i=0; i<mesh.nCells(); i++)
    {
        localIDs[i] = i;
        globalIDs[i] = globalCellMap.toGlobal(i);
    }

    *ierr = ZOLTAN_OK;
}


static int get_mesh_dim(void* data, int* ierr)
{
    const Foam::polyMesh& mesh = *static_cast<const Foam::polyMesh*>(data);

    return mesh.nSolutionD();
}


static void get_geom_list
(
    void* data,
    int nGID,
    int nLID,
    int nCells,
    ZOLTAN_ID_PTR globalIDs,
    ZOLTAN_ID_PTR localIDs,
    int nDim,
    double* cellCentres,
    int* ierr
)
{
    const Foam::polyMesh& mesh = *static_cast<const Foam::polyMesh*>(data);

    if
    (
        (nGID != 1)
     || (nLID != 1)
     || (nCells != mesh.nCells())
     || (nDim != mesh.nSolutionD())
    )
    {
        *ierr = ZOLTAN_FATAL;
        return;
    }

    double* p = cellCentres;

    const Foam::Vector<Foam::label>& sol = mesh.solutionD();
    const Foam::pointField& cc = mesh.cellCentres();

    for (Foam::label celli=0; celli<nCells; celli++)
    {
        const Foam::point& pt = cc[celli];

        for (Foam::direction cmpt=0; cmpt<Foam::vector::nComponents; cmpt++)
        {
            if (sol[cmpt] == 1)
            {
                *p++ = pt[cmpt];
            }
        }
    }

    *ierr = ZOLTAN_OK;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoltanRenumber::zoltanRenumber(const dictionary& renumberDict)
:
    renumberMethod(renumberDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::zoltanRenumber::renumber
(
    const polyMesh& mesh,
    const pointField& points
) const
{
    stringList args(1);
    args[0] = "zoltanRenumber";

    int argc = args.size();
    char* argv[argc];
    for (label i = 0; i < argc; i++)
    {
        argv[i] = strdup(args[i].c_str());
    }

    float ver;
    int rc = Zoltan_Initialize(argc, argv, &ver);

    if (rc != ZOLTAN_OK)
    {
        FatalErrorInFunction
            << "Failed initialising Zoltan" << exit(FatalError);
    }

    struct Zoltan_Struct *zz = Zoltan_Create(PstreamGlobals::MPI_COMM_FOAM);

    {
        // Set default order method to LOCAL_HSFC
        Zoltan_Set_Param(zz, "ORDER_METHOD", "LOCAL_HSFC");
        Zoltan_Set_Param(zz, "ORDER_TYPE", "LOCAL");

        forAllConstIter(IDLList<entry>, coeffsDict_, iter)
        {
            if (!iter().isDict())
            {
                const word& key = iter().keyword();
                const word value(iter().stream());

                Info<< typeName << " : setting parameter " << key
                    << " to " << value << endl;

                Zoltan_Set_Param(zz, key.c_str(), value.c_str());
            }
        }

        void* meshPtr = &const_cast<polyMesh&>(mesh);

        Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, meshPtr);
        Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, meshPtr);

        // Callbacks for geometry
        Zoltan_Set_Num_Geom_Fn(zz, get_mesh_dim, meshPtr);
        Zoltan_Set_Geom_Multi_Fn(zz, get_geom_list, meshPtr);
    }

    // Local to global cell index mapper
    globalIndex globalCellMap(mesh.nCells());

    List<ZOLTAN_ID_TYPE> globalCells(mesh.nCells());
    forAll(globalCells, i)
    {
        globalCells[i] = globalCellMap.toGlobal(i);
    }

    // Old local to new global cell index map
    // Values set by Zoltan_Order
    List<ZOLTAN_ID_TYPE> oldToNew(mesh.nCells());

    // Call Zoltan to reorder the cells
    int err = Zoltan_Order
    (
        zz,
        1,              // int nGID,
        mesh.nCells(),
        globalCells.begin(),
        oldToNew.begin()
    );

    // Check for Zoltan errors
    if (err != ZOLTAN_OK)
    {
        FatalErrorInFunction
            << "Failed Zoltan_Order" << exit(FatalError);
    }

    // Free the argv array
    for (label i = 0; i < argc; i++)
    {
        free(argv[i]);
    }

    // Free the Zoltan_Struct allocated by Zoltan_Create
    Zoltan_Destroy(&zz);

    // Map the oldToNew global cell indices to local
    labelList order(oldToNew.size());
    forAll(order, i)
    {
        order[i] = globalCellMap.toLocal(oldToNew[i]);
    }

    return order;
}


// ************************************************************************* //
