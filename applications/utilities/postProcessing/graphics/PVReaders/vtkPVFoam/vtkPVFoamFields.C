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
#include "vtkPVFoamVolFields.H"
#include "vtkPVFoamPointFields.H"
#include "vtkPVFoamLagrangianFields.H"

// OpenFOAM includes
#include "domainDecomposition.H"
#include "cloud.H"
#include "LagrangianMesh.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::IOobjectList Foam::vtkPVFoam::getObjects
(
    const wordHashSet& selected,
    const fvMesh& mesh,
    const fileName& instance,
    const fileName& local
)
{
    // If nothing is selected then return an empty list
    if (selected.empty())
    {
        return IOobjectList(0);
    }

    // Create the list of objects at the instance
    IOobjectList objects(mesh, instance, local);

    // Add any objects from constant that are not already present
    IOobjectList objectsConstant(mesh, Time::constant(), local);
    forAllIter(IOobjectList, objectsConstant, iter)
    {
        if (!objects.found(iter.key()))
        {
            objects.add
            (
                *objectsConstant.HashPtrTable<IOobject>::remove(iter)
            );
        }
    }

    // Remove everything that is not selected
    forAllIter(IOobjectList, objects, iter)
    {
        if (!selected.found(iter()->name()))
        {
            objects.erase(iter);
        }
    }

    return objects;
}


void Foam::vtkPVFoam::convertFields(vtkMultiBlockDataSet* output)
{
    if (reader_->GetDecomposedCase() && !procDbsPtr_->nProcs()) return;

    const wordHashSet selectedFields =
        getSelected(reader_->GetFieldSelection());

    if (selectedFields.empty()) return;

    // Get objects (fields) for this time - only keep selected fields
    // the region name is already in the mesh db
    IOobjectList objects
    (
        getObjects
        (
            selectedFields,
            reader_->GetDecomposedCase()
          ? procMeshesPtr_->procMeshes().first()
          : procMeshesPtr_->completeMesh(),
            reader_->GetDecomposedCase()
          ? procDbsPtr_->proc0Time().name()
          : procDbsPtr_->completeTime().name()
        )
    );

    if (objects.empty()) return;

    DebugInFunction
        << "Converting OpenFOAM volume fields" << endl;
    forAllConstIter(IOobjectList, objects, iter)
    {
        DebugInfo
            << "    " << iter()->name()
            << " == " << iter()->objectPath(false) << endl;
    }

    PtrList<PrimitivePatchInterpolation<primitivePatch>>
        ppInterpList(procMeshesPtr_->completeMesh().boundaryMesh().size());
    forAll(ppInterpList, i)
    {
        ppInterpList.set
        (
            i,
            new PrimitivePatchInterpolation<primitivePatch>
            (
                procMeshesPtr_->completeMesh().boundaryMesh()[i]
            )
        );
    }

    bool interpFields = reader_->GetInterpolateVolFields();

    #define CONVERT_VOL_FIELDS(Type, nullArg) \
        convertVolFields<Type>                \
        (                                     \
            ppInterpList,                     \
            objects,                          \
            interpFields,                     \
            output                            \
        );
    FOR_ALL_FIELD_TYPES(CONVERT_VOL_FIELDS)
    #undef CONVERT_VOL_FIELDS

    #define CONVERT_VOL_INTERNAL_FIELDS(Type, nullArg) \
        convertVolInternalFields<Type>(objects, output);
    FOR_ALL_FIELD_TYPES(CONVERT_VOL_INTERNAL_FIELDS)
    #undef CONVERT_VOL_INTERNAL_FIELDS

    #define CONVERT_SURFACE_FIELDS(Type, nullArg) \
        convertSurfaceFields<Type>(objects, output);
    FOR_ALL_FIELD_TYPES(CONVERT_SURFACE_FIELDS)
    #undef CONVERT_SURFACE_FIELDS

    #define CONVERT_POINT_FIELDS(Type, nullArg) \
        convertPointFields<Type>(objects, output);
    FOR_ALL_FIELD_TYPES(CONVERT_POINT_FIELDS)
    #undef CONVERT_POINT_FIELDS
}


void Foam::vtkPVFoam::convertlagrangianFields(vtkMultiBlockDataSet* output)
{
    const fileName lagrangianPrefix =
        meshRegion_ == polyMesh::defaultRegion
      ? fileName(lagrangian::cloud::prefix)
      : meshRegion_/lagrangian::cloud::prefix;

    arrayRange& range = arrayRangelagrangian_;

    const wordHashSet selectedFields = getSelected
    (
        reader_->GetlagrangianFieldSelection()
    );

    if (selectedFields.empty()) return;

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word lagrangianName = getPartName(partId);
        const label datasetNo = partDataset_[partId];

        if (!partStatus_[partId] || datasetNo < 0) continue;

        // Find a processor that has some of the cloud in it
        label proci = -1;
        if (reader_->GetDecomposedCase())
        {
            forAll(procDbsPtr_->procTimes(), procj)
            {
                if
                (
                    fileHandler().isDir
                    (
                        fileHandler().filePath
                        (
                            procDbsPtr_->procTimes()[procj].path()
                           /procDbsPtr_().completeTime().name()
                           /lagrangianPrefix
                           /lagrangianName
                        )
                    )
                )
                {
                    proci = procj;
                    break;
                }
            }
        }

        if (reader_->GetDecomposedCase() && proci == -1) continue;

        // Get the lagrangian fields for this time and this cloud
        // but only keep selected fields
        // the region name is already in the mesh db
        IOobjectList objects
        (
            getObjects
            (
                selectedFields,
                reader_->GetDecomposedCase()
              ? procMeshesPtr_->procMeshes()[proci]
              : procMeshesPtr_->completeMesh(),
                procDbsPtr_().completeTime().name(),
                lagrangian::cloud::prefix/lagrangianName
            )
        );

        if (objects.empty()) continue;

        DebugInFunction
            << "Converting OpenFOAM lagrangian fields" << endl;
        forAllConstIter(IOobjectList, objects, iter)
        {
            DebugInfo
                << "    " << iter()->name()
                << " == " << iter()->objectPath(false) << endl;
        }

        #define CONVERT_LAGRANGIAN_FIELDS(Type, nullArg) \
            convertlagrangianFields<Type>(objects, output, datasetNo);
        CONVERT_LAGRANGIAN_FIELDS(label, );
        FOR_ALL_FIELD_TYPES(CONVERT_LAGRANGIAN_FIELDS);
        #undef CONVERT_LAGRANGIAN_FIELDS
    }
}


void Foam::vtkPVFoam::convertLagrangianFields(vtkMultiBlockDataSet* output)
{
    arrayRange& range = arrayRangeLagrangian_;

    const wordHashSet selectedFields = getSelected
    (
        reader_->GetLagrangianFieldSelection()
    );

    if (selectedFields.empty()) return;

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word LagrangianName = getPartName(partId);
        const label datasetNo = partDataset_[partId];

        if (!partStatus_[partId] || datasetNo < 0) continue;

        // Get the Lagrangian fields for this time and this cloud
        // but only keep selected fields
        // the region name is already in the mesh db
        IOobjectList objects
        (
            getObjects
            (
                selectedFields,
                reader_->GetDecomposedCase()
              ? procMeshesPtr_->procMeshes().first()
              : procMeshesPtr_->completeMesh(),
                procDbsPtr_().completeTime().name(),
                LagrangianMesh::prefix/LagrangianName
            )
        );

        if (objects.empty()) continue;

        DebugInFunction
            << "Converting OpenFOAM Lagrangian fields" << endl;
        forAllConstIter(IOobjectList, objects, iter)
        {
            DebugInfo
                << "    " << iter()->name()
                << " == " << iter()->objectPath(false) << endl;
        }

        #define CONVERT_LAGRANGIAN_FIELDS(Type, GeoField) \
            convertLagrangianFields<Type, GeoField>(objects, output, datasetNo);
        CONVERT_LAGRANGIAN_FIELDS(label, LagrangianField);
        FOR_ALL_FIELD_TYPES(CONVERT_LAGRANGIAN_FIELDS, LagrangianField);
        CONVERT_LAGRANGIAN_FIELDS(label, LagrangianInternalField);
        FOR_ALL_FIELD_TYPES(CONVERT_LAGRANGIAN_FIELDS, LagrangianInternalField);
        #undef CONVERT_LAGRANGIAN_FIELDS
    }
}


// ************************************************************************* //
