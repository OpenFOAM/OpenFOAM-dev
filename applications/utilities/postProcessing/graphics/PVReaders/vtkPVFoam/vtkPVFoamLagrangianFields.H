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

InClass
    vtkPVFoam

\*---------------------------------------------------------------------------*/

#ifndef vtkPVFoamLagrangianFields_H
#define vtkPVFoamLagrangianFields_H

#include "vtkPVFoam.H"
#include "IOobjectList.H"
#include "LagrangianFields.H"
#include "vtkOpenFOAMTupleRemap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPVFoam::convertlagrangianFields
(
    const IOobjectList& objects,
    vtkMultiBlockDataSet* output,
    const label datasetNo
)
{
    const arrayRange& range = arrayRangelagrangian_;

    forAllConstIter(IOobjectList, objects, iter)
    {
        // restrict to this IOField<Type>
        if (iter()->headerClassName() == IOField<Type>::typeName)
        {
            IOField<Type> tf(*iter());
            convertlagrangianField(tf, output, range, datasetNo);
        }
    }
}


template<class Type, template<class> class GeoField>
void Foam::vtkPVFoam::convertLagrangianFields
(
    const IOobjectList& objects,
    vtkMultiBlockDataSet* output,
    const label datasetNo,
    const LagrangianMesh& Lmesh
)
{
    const arrayRange& range = arrayRangeLagrangian_;

    forAllConstIter(IOobjectList, objects, iter)
    {
        if (iter()->headerClassName() == GeoField<Type>::typeName)
        {
            GeoField<Type> tf(*iter(), Lmesh);

            convertLagrangianField<Type, GeoField>
            (
                tf,
                output,
                range,
                datasetNo
            );
        }
    }
}


template<class Type>
void Foam::vtkPVFoam::convertlagrangianField
(
    const IOField<Type>& tf,
    vtkMultiBlockDataSet* output,
    const arrayRange& range,
    const label datasetNo
)
{
    const label nComp = pTraits<Type>::nComponents;

    vtkFloatArray* pointData = vtkFloatArray::New();
    pointData->SetNumberOfTuples(tf.size());
    pointData->SetNumberOfComponents(nComp);
    pointData->Allocate(nComp*tf.size());
    pointData->SetName(tf.name().c_str());

    if (debug)
    {
        Info<< "convert lagrangianField: "
            << tf.name()
            << " size = " << tf.size()
            << " nComp=" << nComp
            << " nTuples = " << tf.size() <<  endl;
    }

    float vec[nComp];
    forAll(tf, i)
    {
        const Type& t = tf[i];
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


template<class Type, template<class> class GeoField>
void Foam::vtkPVFoam::convertLagrangianField
(
    const GeoField<Type>& tf,
    vtkMultiBlockDataSet* output,
    const arrayRange& range,
    const label datasetNo
)
{
    static const label nComp = pTraits<Type>::nComponents;

    vtkFloatArray* pointData = vtkFloatArray::New();
    pointData->SetNumberOfTuples(tf.primitiveField().size());
    pointData->SetNumberOfComponents(nComp);
    pointData->Allocate(nComp*tf.primitiveField().size());
    pointData->SetName(tf.name().c_str());

    if (debug)
    {
        Info<< "convert LagrangianField: "
            << tf.name()
            << " size = " << tf.primitiveField().size()
            << " nComp=" << nComp
            << " nTuples = " << tf.primitiveField().size() <<  endl;
    }

    float vec[nComp];

    forAll(tf.primitiveField(), i)
    {
        const Type& t = tf.primitiveField()[i];

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
