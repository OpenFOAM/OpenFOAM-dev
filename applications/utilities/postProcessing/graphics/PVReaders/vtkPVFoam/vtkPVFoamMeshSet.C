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

#include "vtkPVFoam.H"

// OpenFOAM includes
#include "faceSet.H"
#include "pointSet.H"
#include "vtkOpenFOAMPoints.H"

// VTK includes
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vtkPolyData* Foam::vtkPVFoam::faceSetVTKMesh
(
    const fvMesh& mesh,
    const faceSet& fSet
)
{
    vtkPolyData* vtkmesh = vtkPolyData::New();

    if (debug)
    {
        InfoInFunction << endl;
        printMemory();
    }

    // Construct primitivePatch of faces in fSet.

    const faceList& meshFaces = mesh.faces();
    faceList patchFaces(fSet.size());
    label facei = 0;
    forAllConstIter(faceSet, fSet, iter)
    {
        patchFaces[facei++] = meshFaces[iter.key()];
    }
    primitiveFacePatch p(patchFaces, mesh.points());


    // The balance of this routine should be identical to patchVTKMesh

    // Convert OpenFOAM mesh vertices to VTK
    const pointField& points = p.localPoints();

    vtkPoints* vtkpoints = vtkPoints::New();
    vtkpoints->Allocate(points.size());
    forAll(points, i)
    {
        vtkInsertNextOpenFOAMPoint(vtkpoints, points[i]);
    }
    vtkmesh->SetPoints(vtkpoints);
    vtkpoints->Delete();

    // Add faces as polygons
    const faceList& faces = p.localFaces();

    vtkCellArray* vtkcells = vtkCellArray::New();
    vtkcells->Allocate(faces.size());

    forAll(faces, facei)
    {
        const face& f = faces[facei];
        vtkIdType nodeIds[f.size()];

        forAll(f, fp)
        {
            nodeIds[fp] = f[fp];
        }
        vtkcells->InsertNextCell(f.size(), nodeIds);
    }

    vtkmesh->SetPolys(vtkcells);
    vtkcells->Delete();

    if (debug)
    {
        printMemory();
    }

    return vtkmesh;
}


vtkPolyData* Foam::vtkPVFoam::pointSetVTKMesh
(
    const fvMesh& mesh,
    const pointSet& pSet
)
{
    vtkPolyData* vtkmesh = vtkPolyData::New();

    if (debug)
    {
        InfoInFunction << endl;
        printMemory();
    }

    const pointField& meshPoints = mesh.points();

    vtkPoints* vtkpoints = vtkPoints::New();
    vtkpoints->Allocate(pSet.size());

    forAllConstIter(pointSet, pSet, iter)
    {
        vtkInsertNextOpenFOAMPoint(vtkpoints, meshPoints[iter.key()]);
    }

    vtkmesh->SetPoints(vtkpoints);
    vtkpoints->Delete();

    if (debug)
    {
        printMemory();
    }

    return vtkmesh;
}


// ************************************************************************* //
