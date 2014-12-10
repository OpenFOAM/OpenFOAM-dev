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

#ifndef vtkPV3FoamPatchField_H
#define vtkPV3FoamPatchField_H

// VTK includes
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"

#include "vtkOpenFOAMTupleRemap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPV3Foam::convertPatchField
(
    const word& name,
    const Field<Type>& ptf,
    vtkMultiBlockDataSet* output,
    const arrayRange& range,
    const label datasetNo
)
{
    const label nComp = pTraits<Type>::nComponents;

    vtkFloatArray* cellData = vtkFloatArray::New();
    cellData->SetNumberOfTuples(ptf.size());
    cellData->SetNumberOfComponents(nComp);
    cellData->Allocate(nComp*ptf.size());
    cellData->SetName(name.c_str());

    float vec[nComp];
    forAll(ptf, i)
    {
        const Type& t = ptf[i];
        for (direction d=0; d<nComp; ++d)
        {
            vec[d] = component(t, d);
        }
        vtkOpenFOAMTupleRemap<Type>(vec);

        cellData->InsertTuple(i, vec);
    }

    vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, range, datasetNo)
    )   ->GetCellData()
        ->AddArray(cellData);

    cellData->Delete();
}


// as above, but with PointData()
template<class Type>
void Foam::vtkPV3Foam::convertPatchPointField
(
    const word& name,
    const Field<Type>& pptf,
    vtkMultiBlockDataSet* output,
    const arrayRange& range,
    const label datasetNo
)
{
    const label nComp = pTraits<Type>::nComponents;

    vtkFloatArray* pointData = vtkFloatArray::New();
    pointData->SetNumberOfTuples(pptf.size());
    pointData->SetNumberOfComponents(nComp);
    pointData->Allocate(nComp*pptf.size());
    pointData->SetName(name.c_str());

    float vec[nComp];
    forAll(pptf, i)
    {
        const Type& t = pptf[i];
        for (direction d=0; d<nComp; ++d)
        {
            vec[d] = component(t, d);
        }
        vtkOpenFOAMTupleRemap<Type>(vec);

        pointData->InsertTuple(i, vec);
    }

    vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, range, datasetNo)
    )   ->GetPointData()
        ->AddArray(pointData);

    pointData->Delete();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
