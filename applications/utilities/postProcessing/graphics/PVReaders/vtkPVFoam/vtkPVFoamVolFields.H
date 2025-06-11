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

#ifndef vtkPVFoamVolFields_H
#define vtkPVFoamVolFields_H

#include "vtkPVFoam.H"
#include "vtkPVFoamReader.h"
#include "vtkPVFoamSurfaceField.H"
#include "vtkPVFoamPatchField.H"
#include "vtkOpenFOAMTupleRemap.H"

// OpenFOAM includes
#include "domainDecomposition.H"
#include "emptyFvPatchField.H"
#include "faceSet.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "fvFieldReconstructor.H"

// VTK includes
#include "vtkFloatArray.h"
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPVFoam::convertVolFields
(
    const PtrList<PrimitivePatchInterpolation<primitivePatch>>& ppInterpList,
    const IOobjectList& objects,
    const bool interpFields,
    vtkMultiBlockDataSet* output
)
{
    const polyBoundaryMesh& patches =
        procMeshesPtr_->completeMesh().boundaryMesh();

    forAllConstIter(IOobjectList, objects, iter)
    {
        // Restrict to GeometricField<Type, ...>
        if (iter()->headerClassName() != VolField<Type>::typeName) continue;

        // Load the field
        tmp<VolField<Type>> ttf;
        if (reader_->GetDecomposedCase())
        {
            if (!fvReconstructorPtr_.valid())
            {
                fvReconstructorPtr_.set
                (
                    new fvFieldReconstructor
                    (
                        procMeshesPtr_->completeMesh(),
                        procMeshesPtr_->procMeshes(),
                        procMeshesPtr_->procFaceAddressing(),
                        procMeshesPtr_->procCellAddressing(),
                        procMeshesPtr_->procFaceAddressingBf()
                    )
                );
            }

            ttf =
                fvReconstructorPtr_
              ->reconstructVolField<Type>(*iter());
        }
        else
        {
            ttf =
                new VolField<Type>
                (
                    *iter(),
                    procMeshesPtr_->completeMesh()
                );
        }
        const VolField<Type>& tf = ttf();

        // Interpolated field (demand driven)
        autoPtr<PointField<Type>> ptfPtr;
        if (interpFields)
        {
            DebugInFunction
                << "Interpolating field " << tf.name() << endl;

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
            setPolyDecomp_
        );


        // Convert patches - if activated
        for
        (
            int partId = arrayRangePatches_.start();
            partId < arrayRangePatches_.end();
            ++partId
        )
        {
            const word patchName = getPartName(partId);
            const label datasetNo = partDataset_[partId];
            const label patchId = patches.findIndex(patchName);

            if (!partStatus_[partId] || datasetNo < 0 || patchId < 0) continue;

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

        // Convert face zones - if activated
        for
        (
            int partId = arrayRangeFaceZones_.start();
            partId < arrayRangeFaceZones_.end();
            ++partId
        )
        {
            const word zoneName = getPartName(partId);
            const label datasetNo = partDataset_[partId];

            if (!partStatus_[partId] || datasetNo < 0) continue;

            const faceZoneList& zMesh = tf.mesh().faceZones();
            const label zoneId = zMesh.findIndex(zoneName);

            if (zoneId < 0) continue;

            convertSurfaceField
            (
                tf,
                output,
                arrayRangeFaceZones_,
                datasetNo,
                tf.mesh(),
                zMesh[zoneId]
            );
        }

        // Convert face sets - if activated
        for
        (
            int partId = arrayRangeFaceSets_.start();
            partId < arrayRangeFaceSets_.end();
            ++partId
        )
        {
            const word selectName = getPartName(partId);
            const label datasetNo = partDataset_[partId];

            if (!partStatus_[partId] || datasetNo < 0) continue;

            const autoPtr<faceSet> fSetPtr =
                reader_->GetDecomposedCase()
              ? procMeshesPtr_->reconstructSet<faceSet>(selectName)
              : autoPtr<faceSet>(new faceSet(tf.mesh(), selectName));

            convertSurfaceField
            (
                tf,
                output,
                arrayRangeFaceSets_,
                datasetNo,
                tf.mesh(),
                fSetPtr().toc()
            );
        }
    }
}


template<class Type>
void Foam::vtkPVFoam::convertVolInternalFields
(
    const IOobjectList& objects,
    vtkMultiBlockDataSet* output
)
{
    forAllConstIter(IOobjectList, objects, iter)
    {
        // Restrict to GeometricField<Type, ...>::Internal
        if (iter()->headerClassName() != VolInternalField<Type>::typeName)
            continue;

        // Load the field
        tmp<VolInternalField<Type>> ttf;
        if (reader_->GetDecomposedCase())
        {
            if (!fvReconstructorPtr_.valid())
            {
                fvReconstructorPtr_.set
                (
                    new fvFieldReconstructor
                    (
                        procMeshesPtr_->completeMesh(),
                        procMeshesPtr_->procMeshes(),
                        procMeshesPtr_->procFaceAddressing(),
                        procMeshesPtr_->procCellAddressing(),
                        procMeshesPtr_->procFaceAddressingBf()
                    )
                );
            }

            ttf =
                fvReconstructorPtr_
              ->reconstructVolInternalField<Type>(*iter());
        }
        else
        {
            ttf =
                new VolInternalField<Type>
                (
                    *iter(),
                    procMeshesPtr_->completeMesh()
                );
        }
        const VolInternalField<Type>& tf = ttf();

        // Convert activated internalMesh regions
        convertVolInternalFieldBlock<Type>
        (
            tf,
            output,
            arrayRangeVolume_,
            regionPolyDecomp_
        );

        // Convert activated cellZones
        convertVolInternalFieldBlock<Type>
        (
            tf,
            output,
            arrayRangeCellZones_,
            zonePolyDecomp_
        );

        // Convert activated cellSets
        convertVolInternalFieldBlock<Type>
        (
            tf,
            output,
            arrayRangeCellSets_,
            setPolyDecomp_
        );
    }
}


template<class Type>
void Foam::vtkPVFoam::convertVolFieldBlock
(
    const VolField<Type>& tf,
    autoPtr<PointField<Type>>& ptfPtr,
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
            convertVolInternalField<Type>
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
void Foam::vtkPVFoam::convertVolInternalFieldBlock
(
    const VolInternalField<Type>& tf,
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
            convertVolInternalField<Type>
            (
                tf,
                output,
                range,
                datasetNo,
                decompLst[datasetNo]
            );
        }
    }
}


template<class Type>
void Foam::vtkPVFoam::convertVolInternalField
(
    const VolInternalField<Type>& tf,
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

    DebugInFunction
        << "Converting Vol field: " << tf.name()
        << " size=" << tf.size() << " (" << superCells.size()
        << "), nComp=" << nComp << endl;

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
