/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "vtkPVFoam.H"

// OpenFOAM includes
#include "Cloud.H"
#include "LagrangianMesh.H"
#include "fvMesh.H"
#include "IOobjectList.H"
#include "passiveParticle.H"
#include "vtkOpenFOAMPoints.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vtkPolyData* Foam::vtkPVFoam::lagrangianVTKMesh
(
    const fvMesh& mesh,
    const word& cloudName
)
{
    vtkPolyData* vtkmesh = nullptr;

    if (debug)
    {
        InfoInFunction
            << "timePath "
            << mesh.time().timePath()/lagrangian::cloud::prefix/cloudName
            << endl;
        printMemory();
    }

    // The region name is already in the mesh db
    IOobjectList cloudObjs
    (
        mesh,
        mesh.time().name(),
        lagrangian::cloud::prefix/cloudName
    );

    IOobject* positionsPtr = cloudObjs.lookup(word("positions"));
    if (positionsPtr)
    {
        lagrangian::Cloud<passiveParticle> parcels(mesh, cloudName, false);

        if (debug)
        {
            Info<< "    cloud with " << parcels.size()
                << " parcels" << endl;
        }

        vtkmesh = vtkPolyData::New();
        vtkPoints* vtkpoints = vtkPoints::New();
        vtkCellArray* vtkcells = vtkCellArray::New();

        vtkpoints->Allocate(parcels.size());
        vtkcells->Allocate(parcels.size());

        vtkIdType particleId = 0;
        forAllConstIter(lagrangian::Cloud<passiveParticle>, parcels, iter)
        {
            vtkInsertNextOpenFOAMPoint(vtkpoints, iter().position(mesh));

            vtkcells->InsertNextCell(1, &particleId);
            particleId++;
        }

        vtkmesh->SetPoints(vtkpoints);
        vtkpoints->Delete();

        vtkmesh->SetVerts(vtkcells);
        vtkcells->Delete();
    }

    if (debug)
    {
        printMemory();
    }

    return vtkmesh;
}


vtkPolyData* Foam::vtkPVFoam::LagrangianVTKMesh
(
    const fvMesh& mesh,
    const word& LagrangianName,
    autoPtr<LagrangianMesh>& LmeshPtr
)
{
    vtkPolyData* vtkmesh = nullptr;

    if (debug)
    {
        InfoInFunction
            << "timePath "
            << mesh.time().timePath()/LagrangianMesh::prefix/LagrangianName
            << endl;
        printMemory();
    }

    LmeshPtr.reset(new LagrangianMesh(mesh, LagrangianName));

    const LagrangianVectorInternalField Lpositions(LmeshPtr->position());

    vtkmesh = vtkPolyData::New();

    vtkPoints* vtkpoints = vtkPoints::New();
    vtkCellArray* vtkcells = vtkCellArray::New();

    vtkpoints->Allocate(LmeshPtr->size());
    vtkcells->Allocate(LmeshPtr->size());

    for (vtkIdType i = 0; i < LmeshPtr->size(); ++ i)
    {
        vtkInsertNextOpenFOAMPoint(vtkpoints, Lpositions[i]);

        vtkcells->InsertNextCell(1, &i);
    }

    vtkmesh->SetPoints(vtkpoints);
    vtkpoints->Delete();

    vtkmesh->SetVerts(vtkcells);
    vtkcells->Delete();

    if (debug)
    {
        printMemory();
    }

    return vtkmesh;
}


// ************************************************************************* //
