/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

InClass
    vtkPV3Foam

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3FoamFaceField_H
#define vtkPV3FoamFaceField_H

// VTK includes
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"

#include "vtkOpenFOAMTupleRemap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPV3Foam::convertFaceField
(
    const GeometricField<Type, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output,
    const arrayRange& range,
    const label datasetNo,
    const fvMesh& mesh,
    const labelList& faceLabels
)
{
    const label nComp = pTraits<Type>::nComponents;
    const label nInternalFaces = mesh.nInternalFaces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeigh = mesh.faceNeighbour();

    vtkFloatArray* cellData = vtkFloatArray::New();
    cellData->SetNumberOfTuples(faceLabels.size());
    cellData->SetNumberOfComponents(nComp);
    cellData->Allocate(nComp*faceLabels.size());
    cellData->SetName(tf.name().c_str());

    if (debug)
    {
        Info<< "convert convertFaceField: "
            << tf.name()
            << " size = " << tf.size()
            << " nComp=" << nComp
            << " nTuples = " << faceLabels.size() <<  endl;
    }

    float vec[nComp];

    // for interior faces: average owner/neighbour
    // for boundary faces: owner
    forAll(faceLabels, faceI)
    {
        const label faceNo = faceLabels[faceI];
        if (faceNo < nInternalFaces)
        {
            Type t = 0.5*(tf[faceOwner[faceNo]] + tf[faceNeigh[faceNo]]);

            for (direction d=0; d<nComp; ++d)
            {
                vec[d] = component(t, d);
            }
        }
        else
        {
            const Type& t = tf[faceOwner[faceNo]];
            for (direction d=0; d<nComp; ++d)
            {
                vec[d] = component(t, d);
            }
        }
        vtkOpenFOAMTupleRemap<Type>(vec);

        cellData->InsertTuple(faceI, vec);
    }


    vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, range, datasetNo)
    )   ->GetCellData()
        ->AddArray(cellData);

    cellData->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
