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
#include "vtkPVFoamReader.h"
#include "vtkOpenFOAMPoints.H"

// OpenFOAM includes
#include "domainDecomposition.H"
#include "lagrangianFieldReconstructor.H"
#include "LagrangianMesh.H"
#include "LagrangianFieldReconstructor.H"
#include "passiveParticleCloud.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vtkPolyData* Foam::vtkPVFoam::lagrangianVTKMesh
(
    const word& lagrangianName,
    autoPtr<lagrangianFieldReconstructor>& lreconstructorPtr
)
{
    vtkPolyData* vtkmesh = nullptr;

    if (reader_->GetDecomposedCase())
    {
        lreconstructorPtr.reset
        (
            new lagrangianFieldReconstructor
            (
                procMeshesPtr_->completeMesh(),
                procMeshesPtr_->procMeshes(),
                procMeshesPtr_->procFaceAddressing(),
                procMeshesPtr_->procCellAddressing(),
                lagrangianName
            )
        );
    }

    autoPtr<passiveParticleCloud> tcloud
    (
        reader_->GetDecomposedCase()
      ? lreconstructorPtr->completeCloud().ptr()
      : new passiveParticleCloud
        (
            procMeshesPtr_->completeMesh(),
            lagrangianName,
            false
        )
    );
    const passiveParticleCloud& cloud = tcloud();

    vtkmesh = vtkPolyData::New();
    vtkPoints* vtkpoints = vtkPoints::New();
    vtkCellArray* vtkcells = vtkCellArray::New();

    vtkpoints->Allocate(cloud.size());
    vtkcells->Allocate(cloud.size());

    vtkIdType particleId = 0;
    forAllConstIter(lagrangian::Cloud<passiveParticle>, cloud, iter)
    {
        vtkInsertNextOpenFOAMPoint(vtkpoints, iter().position(cloud.pMesh()));

        vtkcells->InsertNextCell(1, &particleId);
        particleId++;
    }

    vtkmesh->SetPoints(vtkpoints);
    vtkpoints->Delete();

    vtkmesh->SetVerts(vtkcells);
    vtkcells->Delete();

    return vtkmesh;
}


vtkPolyData* Foam::vtkPVFoam::LagrangianVTKMesh
(
    const word& LagrangianName,
    autoPtr<LagrangianMesh>& LmeshPtr,
    autoPtr<LagrangianFieldReconstructor>& LreconstructorPtr
)
{
    vtkPolyData* vtkmesh = nullptr;

    if (reader_->GetDecomposedCase())
    {
        LreconstructorPtr.reset
        (
            new LagrangianFieldReconstructor
            (
                procMeshesPtr_->completeMesh(),
                procMeshesPtr_->procMeshes(),
                procMeshesPtr_->procFaceAddressing(),
                procMeshesPtr_->procCellAddressing(),
                LagrangianName
            )
        );
    }
    else
    {
        LmeshPtr.reset
        (
            new LagrangianMesh(procMeshesPtr_->completeMesh(), LagrangianName)
        );
    }

    const LagrangianMesh& Lmesh =
        reader_->GetDecomposedCase()
      ? LreconstructorPtr().completeMesh()
      : LmeshPtr();

    const LagrangianInternalVectorField Lpositions(Lmesh.position());

    vtkmesh = vtkPolyData::New();

    vtkPoints* vtkpoints = vtkPoints::New();
    vtkCellArray* vtkcells = vtkCellArray::New();

    vtkpoints->Allocate(Lmesh.size());
    vtkcells->Allocate(Lmesh.size());

    for (vtkIdType i = 0; i < Lmesh.size(); ++ i)
    {
        vtkInsertNextOpenFOAMPoint(vtkpoints, Lpositions[i]);

        vtkcells->InsertNextCell(1, &i);
    }

    vtkmesh->SetPoints(vtkpoints);
    vtkpoints->Delete();

    vtkmesh->SetVerts(vtkcells);
    vtkcells->Delete();

    return vtkmesh;
}


// ************************************************************************* //
