/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

InClass
    vtkPVFoam

\*---------------------------------------------------------------------------*/

#ifndef vtkPVFoamPointFields_H
#define vtkPVFoamPointFields_H

// OpenFOAM includes
#include "interpolatePointToCell.H"

#include "vtkOpenFOAMTupleRemap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPVFoam::convertPointFields
(
    const fvMesh& mesh,
    const pointMesh& pMesh,
    const IOobjectList& objects,
    vtkMultiBlockDataSet* output
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAllConstIter(IOobjectList, objects, iter)
    {
        const word& fieldName = iter()->name();
        // restrict to this GeometricField<Type, ...>
        if
        (
            iter()->headerClassName()
         != GeometricField<Type, pointPatchField, pointMesh>::typeName
        )
        {
            continue;
        }

        if (debug)
        {
            Info<< "Foam::vtkPVFoam::convertPointFields : "
                << fieldName << endl;
        }

        GeometricField<Type, pointPatchField, pointMesh> ptf
        (
            *iter(),
            pMesh
        );


        // Convert activated internalMesh regions
        convertPointFieldBlock
        (
            ptf,
            output,
            arrayRangeVolume_,
            regionPolyDecomp_
        );

        // Convert activated cellZones
        convertPointFieldBlock
        (
            ptf,
            output,
            arrayRangeCellZones_,
            zonePolyDecomp_
        );

        // Convert activated cellSets
        convertPointFieldBlock
        (
            ptf,
            output,
            arrayRangeCellSets_,
            csetPolyDecomp_
        );


        //
        // Convert patches - if activated
        //
        for
        (
            int partId = arrayRangePatches_.start();
            partId < arrayRangePatches_.end();
            ++partId
        )
        {
            const word patchName = getPartName(partId);
            const label datasetNo = partDataset_[partId];
            const label patchId = patches.findPatchID(patchName);

            if (!partStatus_[partId] || datasetNo < 0 || patchId < 0)
            {
                continue;
            }

            convertPatchPointField
            (
                fieldName,
                ptf.boundaryField()[patchId].patchInternalField()(),
                output,
                arrayRangePatches_,
                datasetNo
            );
        }

        //
        // Convert faceZones - if activated
        //
        for
        (
            int partId = arrayRangeFaceZones_.start();
            partId < arrayRangeFaceZones_.end();
            ++partId
        )
        {
            const word zoneName = getPartName(partId);
            const label datasetNo = partDataset_[partId];
            const label zoneId = mesh.faceZones().findZoneID(zoneName);

            if (!partStatus_[partId] || datasetNo < 0 || zoneId < 0)
            {
                continue;
            }

            // Extract the field on the zone
            Field<Type> fld
            (
                ptf.primitiveField(),
                mesh.faceZones()[zoneId]().meshPoints()
            );

            convertPatchPointField
            (
                fieldName,
                fld,
                output,
                arrayRangeFaceZones_,
                datasetNo
            );
        }
    }
}


template<class Type>
void Foam::vtkPVFoam::convertPointFieldBlock
(
    const GeometricField<Type, pointPatchField, pointMesh>& ptf,
    vtkMultiBlockDataSet* output,
    const arrayRange& range,
    const List<polyDecomp>& decompLst
)
{
   for (int partId = range.start(); partId < range.end(); ++partId)
   {
       const label datasetNo = partDataset_[partId];

       if (datasetNo >= 0 && partStatus_[partId])
       {
           convertPointField
           (
               ptf,
               GeometricField<Type, fvPatchField, volMesh>::null(),
               output,
               range,
               datasetNo,
               decompLst[datasetNo]
           );
       }
   }
}


template<class Type>
void Foam::vtkPVFoam::convertPointField
(
    const GeometricField<Type, pointPatchField, pointMesh>& ptf,
    const GeometricField<Type, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output,
    const arrayRange& range,
    const label datasetNo,
    const polyDecomp& decomp
)
{
    const label nComp = pTraits<Type>::nComponents;
    const labelList& addPointCellLabels = decomp.addPointCellLabels();
    const labelList& pointMap = decomp.pointMap();

    // use a pointMap or address directly into mesh
    label nPoints;
    if (pointMap.size())
    {
        nPoints = pointMap.size();
    }
    else
    {
        nPoints = ptf.size();
    }

    vtkFloatArray* pointData = vtkFloatArray::New();
    pointData->SetNumberOfTuples(nPoints + addPointCellLabels.size());
    pointData->SetNumberOfComponents(nComp);
    pointData->Allocate(nComp*(nPoints + addPointCellLabels.size()));

    // Note: using the name of the original volField
    // not the name generated by the interpolation "volPointInterpolate(<name>)"

    if (&tf != &GeometricField<Type, fvPatchField, volMesh>::null())
    {
        pointData->SetName(tf.name().c_str());
    }
    else
    {
        pointData->SetName(ptf.name().c_str());
    }

    if (debug)
    {
        Info<< "convert convertPointField: "
            << ptf.name()
            << " size = " << nPoints
            << " nComp=" << nComp
            << " nTuples = " << (nPoints + addPointCellLabels.size())
            <<  endl;
    }

    float vec[nComp];

    if (pointMap.size())
    {
        forAll(pointMap, i)
        {
            const Type& t = ptf[pointMap[i]];
            for (direction d=0; d<nComp; ++d)
            {
                vec[d] = component(t, d);
            }
            vtkOpenFOAMTupleRemap<Type>(vec);

            pointData->InsertTuple(i, vec);
        }
    }
    else
    {
        forAll(ptf, i)
        {
            const Type& t = ptf[i];
            for (direction d=0; d<nComp; ++d)
            {
                vec[d] = component(t, d);
            }
            vtkOpenFOAMTupleRemap<Type>(vec);

            pointData->InsertTuple(i, vec);
        }
    }

    // continue insertion from here
    label i = nPoints;

    if (&tf != &GeometricField<Type, fvPatchField, volMesh>::null())
    {
        forAll(addPointCellLabels, apI)
        {
            const Type& t = tf[addPointCellLabels[apI]];
            for (direction d=0; d<nComp; ++d)
            {
                vec[d] = component(t, d);
            }
            vtkOpenFOAMTupleRemap<Type>(vec);

            pointData->InsertTuple(i++, vec);
        }
    }
    else
    {
        forAll(addPointCellLabels, apI)
        {
            Type t = interpolatePointToCell(ptf, addPointCellLabels[apI]);
            for (direction d=0; d<nComp; ++d)
            {
                vec[d] = component(t, d);
            }
            vtkOpenFOAMTupleRemap<Type>(vec);

            pointData->InsertTuple(i++, vec);
        }
    }

    vtkUnstructuredGrid::SafeDownCast
    (
        GetDataSetFromBlock(output, range, datasetNo)
    )   ->GetPointData()
        ->AddArray(pointData);

    pointData->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
