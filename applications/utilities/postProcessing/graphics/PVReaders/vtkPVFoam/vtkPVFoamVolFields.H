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

#ifndef vtkPVFoamVolFields_H
#define vtkPVFoamVolFields_H

// OpenFOAM includes
#include "emptyFvPatchField.H"
#include "wallPolyPatch.H"
#include "faceSet.H"
#include "volPointInterpolation.H"

#include "vtkPVFoamFaceField.H"
#include "vtkPVFoamPatchField.H"

#include "vtkOpenFOAMTupleRemap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPVFoam::convertVolFields
(
    const fvMesh& mesh,
    const PtrList<PrimitivePatchInterpolation<primitivePatch>>& ppInterpList,
    const IOobjectList& objects,
    const bool interpFields,
    vtkMultiBlockDataSet* output
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAllConstIter(IOobjectList, objects, iter)
    {
        // restrict to GeometricField<Type, ...>
        if
        (
            iter()->headerClassName()
         != GeometricField<Type, fvPatchField, volMesh>::typeName
        )
        {
            continue;
        }

        // Load field
        GeometricField<Type, fvPatchField, volMesh> tf
        (
            *iter(),
            mesh
        );

        // Interpolated field (demand driven)
        autoPtr<GeometricField<Type, pointPatchField, pointMesh>> ptfPtr;
        if (interpFields)
        {
            if (debug)
            {
                Info<< "convertVolFieldBlock interpolating:" << tf.name()
                    << endl;
            }

            ptfPtr.reset
            (
                volPointInterpolation::New(tf.mesh()).interpolate(tf).ptr()
            );
        }


        // Convert activated internalMesh regions
        convertVolFieldBlock
        (
            tf,
            ptfPtr,
            output,
            arrayRangeVolume_,
            regionPolyDecomp_
        );

        // Convert activated cellZones
        convertVolFieldBlock
        (
            tf,
            ptfPtr,
            output,
            arrayRangeCellZones_,
            zonePolyDecomp_
        );

        // Convert activated cellSets
        convertVolFieldBlock
        (
            tf,
            ptfPtr,
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

            const fvPatchField<Type>& ptf = tf.boundaryField()[patchId];

            if
            (
                isType<emptyFvPatchField<Type>>(ptf)
             ||
                (
                    reader_->GetExtrapolatePatches()
                && !polyPatch::constraintType(patches[patchId].type())
                )
            )
            {
                fvPatch p(ptf.patch().patch(), tf.mesh().boundary());

                tmp<Field<Type>> tpptf
                (
                    fvPatchField<Type>(p, tf).patchInternalField()
                );

                convertPatchField
                (
                    tf.name(),
                    tpptf(),
                    output,
                    arrayRangePatches_,
                    datasetNo
                );

                if (interpFields)
                {
                    convertPatchPointField
                    (
                        tf.name(),
                        ppInterpList[patchId].faceToPointInterpolate(tpptf)(),
                        output,
                        arrayRangePatches_,
                        datasetNo
                    );
                }
            }
            else
            {
                convertPatchField
                (
                    tf.name(),
                    ptf,
                    output,
                    arrayRangePatches_,
                    datasetNo
                );

                if (interpFields)
                {
                    convertPatchPointField
                    (
                        tf.name(),
                        ppInterpList[patchId].faceToPointInterpolate(ptf)(),
                        output,
                        arrayRangePatches_,
                        datasetNo
                    );
                }
            }
        }

        //
        // Convert face zones - if activated
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

            if (!partStatus_[partId] || datasetNo < 0)
            {
                continue;
            }

            const faceZoneMesh& zMesh = mesh.faceZones();
            const label zoneId = zMesh.findZoneID(zoneName);

            if (zoneId < 0)
            {
                continue;
            }

            convertFaceField
            (
                tf,
                output,
                arrayRangeFaceZones_,
                datasetNo,
                mesh,
                zMesh[zoneId]
            );

            // TODO: points
        }

        //
        // Convert face sets - if activated
        //
        for
        (
            int partId = arrayRangeFaceSets_.start();
            partId < arrayRangeFaceSets_.end();
            ++partId
        )
        {
            const word selectName = getPartName(partId);
            const label datasetNo = partDataset_[partId];

            if (!partStatus_[partId] || datasetNo < 0)
            {
                continue;
            }

            const faceSet fSet(mesh, selectName);

            convertFaceField
            (
                tf,
                output,
                arrayRangeFaceSets_,
                datasetNo,
                mesh,
                fSet.toc()
            );

            // TODO: points
        }
    }
}


template<class Type>
void Foam::vtkPVFoam::convertVolFieldBlock
(
    const GeometricField<Type, fvPatchField, volMesh>& tf,
    autoPtr<GeometricField<Type, pointPatchField, pointMesh>>& ptfPtr,
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
            convertVolField
            (
                tf,
                output,
                range,
                datasetNo,
                decompLst[datasetNo]
            );

            if (ptfPtr.valid())
            {
                convertPointField
                (
                    ptfPtr(),
                    tf,
                    output,
                    range,
                    datasetNo,
                    decompLst[datasetNo]
                );
            }
        }
    }
}


template<class Type>
void Foam::vtkPVFoam::convertVolField
(
    const GeometricField<Type, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output,
    const arrayRange& range,
    const label datasetNo,
    const polyDecomp& decompInfo
)
{
    const label nComp = pTraits<Type>::nComponents;
    const labelList& superCells = decompInfo.superCells();

    vtkFloatArray* celldata = vtkFloatArray::New();
    celldata->SetNumberOfTuples(superCells.size());
    celldata->SetNumberOfComponents(nComp);
    celldata->Allocate(nComp*superCells.size());
    celldata->SetName(tf.name().c_str());

    if (debug)
    {
        Info<< "convert volField: "
            << tf.name()
            << " size = " << tf.size()
            << " nComp=" << nComp
            << " nTuples = " << superCells.size() <<  endl;
    }

    float vec[nComp];
    forAll(superCells, i)
    {
        const Type& t = tf[superCells[i]];
        for (direction d=0; d<nComp; ++d)
        {
            vec[d] = component(t, d);
        }
        vtkOpenFOAMTupleRemap<Type>(vec);

        celldata->InsertTuple(i, vec);
    }

    vtkUnstructuredGrid::SafeDownCast
    (
        GetDataSetFromBlock(output, range, datasetNo)
    )   ->GetCellData()
        ->AddArray(celldata);

    celldata->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
