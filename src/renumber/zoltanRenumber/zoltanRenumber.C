/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Class
    Foam::zoltanRenumber

Description
    Renumber using Zoltan

    Zoltan install:
    - in your ~/.bashrc:
            export ZOLTAN_ARCH_DIR=\
                $WM_THIRD_PARTY_DIR/platforms/linux64Gcc/Zoltan_XXX
    - unpack into $WM_THIRD_PARTY_DIR
    - cd Zoltan_XXX
    - mkdir build
    - cd build
    - export CCFLAGS="-fPIC"
    - export CXXFLAGS="-fPIC"
    - export CFLAGS="-fPIC"
    - export LDFLAGS="-shared"
    - ../configure \
        --prefix=$ZOLTAN_ARCH_DIR \
        --with-ccflags=-fPIC --with-cxxflags=-fPIC --with-ldflags=-shared

SourceFiles
    zoltanRenumber.C

\*---------------------------------------------------------------------------*/

#include "zoltanRenumber.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "labelIOList.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "globalIndex.H"
#include "uint.H"

#include "zoltan.h"
#include <mpi.h>

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


static void get_vertex_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr)
{
    const Foam::polyMesh& mesh = *static_cast<const Foam::polyMesh*>(data);

    *ierr = ZOLTAN_OK;

  /* In this example, return the IDs of our vertices, but no weights.
   * Zoltan will assume equally weighted vertices.
   */

    wgt_dim = 0;
    obj_wgts = nullptr;

    for (Foam::label i=0; i<mesh.nCells(); i++)
    {
        globalID[i] = i;        // should be global
        localID[i] = i;
    }
}


static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int *numEdges, int *ierr)
{
  const Foam::polyMesh& mesh = *static_cast<const Foam::polyMesh*>(data);

    if ((sizeGID != 1) || (sizeLID != 1) || (num_obj != mesh.nCells()))
    {
        *ierr = ZOLTAN_FATAL;
        return;
    }

    for (Foam::label i=0; i < num_obj ;i++)
    {
        Foam::label celli = localID[i];
        const Foam::cell& cFaces = mesh.cells()[celli];
        forAll(cFaces, cFacei)
        {
            Foam::label n = 0;
            if (mesh.isInternalFace(cFaces[cFacei]))
            {
                n++;
            }
            numEdges[i] = n;
        }
    }

    *ierr = ZOLTAN_OK;
}


static void get_edge_list(void *data, int sizeGID, int sizeLID,
        int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int *num_edges,
        ZOLTAN_ID_PTR nborGID, int *nborProc,
        int wgt_dim, float *ewgts, int *ierr)
{
    const Foam::polyMesh& mesh = *static_cast<const Foam::polyMesh*>(data);

    if
    (
        (sizeGID != 1)
     || (sizeLID != 1)
     || (num_obj != mesh.nCells())
     || (wgt_dim != 1)
    )
    {
        *ierr = ZOLTAN_FATAL;
        return;
    }

    ZOLTAN_ID_TYPE* nextNbor = nborGID;
    int* nextProc = nborProc;
    float* nextWgt = ewgts;

    for (Foam::label i=0; i < num_obj; i++)
    {
        Foam::label celli = localID[i];

        const Foam::cell& cFaces = mesh.cells()[celli];
        forAll(cFaces, cFacei)
        {
            Foam::label n = 0;

            Foam::label facei = cFaces[cFacei];
            if (mesh.isInternalFace(facei))
            {
                Foam::label nbr = mesh.faceOwner()[facei];
                if (nbr == celli)
                {
                    nbr = mesh.faceNeighbour()[facei];
                }

                // Note: global index
                *nextNbor++ = nbr;
                *nextProc++ = 0;
                *nextWgt++ = 1.0;

                n++;
            }
            if (n != num_edges[i])
            {
                *ierr = ZOLTAN_FATAL;
                return;
            }
        }
    }
    *ierr = ZOLTAN_OK;
}


static int get_mesh_dim(void *data, int *ierr)
{
    const Foam::polyMesh& mesh = *static_cast<const Foam::polyMesh*>(data);

    return mesh.nSolutionD();
}


static void get_geom_list
(
    void *data,
    int num_gid_entries,
    int num_lid_entries,
    int num_obj,
    ZOLTAN_ID_PTR global_ids,
    ZOLTAN_ID_PTR local_ids,
    int num_dim,
    double *geom_vec,
    int *ierr
)
{
    const Foam::polyMesh& mesh = *static_cast<const Foam::polyMesh*>(data);

    if
    (
        (num_gid_entries != 1)
     || (num_lid_entries != 1)
     || (num_obj != mesh.nCells())
     || (num_dim != mesh.nSolutionD())
    )
    {
        *ierr = ZOLTAN_FATAL;
        return;
    }

    double* p = geom_vec;


    const Foam::Vector<Foam::label>& sol = mesh.solutionD();

    const Foam::pointField& cc = mesh.cellCentres();

    for (Foam::label celli = 0; celli < num_obj; celli++)
    {
        const Foam::point& pt = cc[celli];

        for (Foam::direction cmpt = 0; cmpt < Foam::vector::nComponents; cmpt++)
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
    renumberMethod(renumberDict),
    coeffsDict_(renumberDict.optionalSubDict(typeName+"Coeffs"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::zoltanRenumber::renumber
(
    const polyMesh& pMesh,
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

    Foam::Pout<< "Initialised to " << ver << Foam::endl;

    if (rc != ZOLTAN_OK)
    {
        FatalErrorInFunction
            << "Failed initialising Zoltan" << exit(FatalError);
    }

    struct Zoltan_Struct *zz = Zoltan_Create(PstreamGlobals::MPI_COMM_FOAM);

    polyMesh& mesh = const_cast<polyMesh&>(pMesh);


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


    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, &mesh);
    Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, &mesh);

    // Callbacks for geometry
    Zoltan_Set_Num_Geom_Fn(zz, get_mesh_dim, &mesh);
    Zoltan_Set_Geom_Multi_Fn(zz, get_geom_list, &mesh);

    // Callbacks for connectivity
    Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list, &mesh);
    Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list, &mesh);



    //Note: !global indices
    List<ZOLTAN_ID_TYPE> wantedCells(mesh.nCells());

    globalIndex globalCells(mesh.nCells());
    forAll(wantedCells, i)
    {
        //wantedCells[i] = i;
        wantedCells[i] = globalCells.toGlobal(i);
    }

    List<ZOLTAN_ID_TYPE> oldToNew(mesh.nCells());

    int err = Zoltan_Order
    (
        zz,
        1,                                //int num_gid_entries,
        mesh.globalData().nTotalCells(),  //int num_obj,
        wantedCells.begin(),
        oldToNew.begin()
    );

    if (err != ZOLTAN_OK)
    {
        FatalErrorInFunction
            << "Failed Zoltan_Order" << exit(FatalError);
    }


    for (label i = 0; i < argc; i++)
    {
        free(argv[i]);
    }


    labelList order(oldToNew.size());
    forAll(order, i)
    {
        order[i] = oldToNew[i];
    }
    return order;
}


// ************************************************************************* //
