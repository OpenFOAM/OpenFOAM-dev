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

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"

// OpenFOAM includes
#include "cloud.H"
#include "LagrangianMesh.H"
#include "IOobjectList.H"
#include "vtkPVFoam.H"
#include "vtkPVFoamReader.h"
#include "vtkPVFoamVolFields.H"
#include "vtkPVFoamPointFields.H"
#include "vtkPVFoamLagrangianFields.H"


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
    IOobjectList objectsConstant(mesh, dbPtr_().constant(), local);
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
    const fvMesh& mesh = *meshPtr_;

    wordHashSet selectedFields = getSelected
    (
        reader_->GetFieldSelection()
    );

    if (selectedFields.empty())
    {
        return;
    }

    // Get objects (fields) for this time - only keep selected fields
    // the region name is already in the mesh db
    IOobjectList objects
    (
        getObjects
        (
            selectedFields,
            mesh,
            dbPtr_().name()
        )
    );

    if (objects.empty())
    {
        return;
    }

    if (debug)
    {
        InfoInFunction<< nl
            << "    converting OpenFOAM volume fields" << endl;
        forAllConstIter(IOobjectList, objects, iter)
        {
            Info<< "  " << iter()->name()
                << " == " << iter()->objectPath(false) << nl;
        }
        printMemory();
    }


    PtrList<PrimitivePatchInterpolation<primitivePatch>>
        ppInterpList(mesh.boundaryMesh().size());

    forAll(ppInterpList, i)
    {
        ppInterpList.set
        (
            i,
            new PrimitivePatchInterpolation<primitivePatch>
            (
                mesh.boundaryMesh()[i]
            )
        );
    }


    bool interpFields = reader_->GetInterpolateVolFields();

    convertVolFields<scalar>
    (
        mesh, ppInterpList, objects, interpFields, output
    );
    convertVolFields<vector>
    (
        mesh, ppInterpList, objects, interpFields, output
    );
    convertVolFields<sphericalTensor>
    (
        mesh, ppInterpList, objects, interpFields, output
    );
    convertVolFields<symmTensor>
    (
        mesh, ppInterpList, objects, interpFields, output
    );
    convertVolFields<tensor>
    (
        mesh, ppInterpList, objects, interpFields, output
    );

    convertVolInternalFields<scalar>
    (
        mesh, objects, output
    );
    convertVolInternalFields<vector>
    (
        mesh, objects, output
    );
    convertVolInternalFields<sphericalTensor>
    (
        mesh, objects, output
    );
    convertVolInternalFields<symmTensor>
    (
        mesh, objects, output
    );
    convertVolInternalFields<tensor>
    (
        mesh, objects, output
    );

    convertSurfaceFields<scalar>(mesh, objects, output);
    convertSurfaceFields<vector>(mesh, objects, output);
    convertSurfaceFields<sphericalTensor>(mesh, objects, output);
    convertSurfaceFields<symmTensor>(mesh, objects, output);
    convertSurfaceFields<tensor>(mesh, objects, output);

    // Construct interpolation on the raw mesh
    const pointMesh& pMesh = pointMesh::New(mesh);

    convertPointFields<scalar>(mesh, pMesh, objects, output);
    convertPointFields<vector>(mesh, pMesh, objects, output);
    convertPointFields<sphericalTensor>(mesh, pMesh, objects, output);
    convertPointFields<symmTensor>(mesh, pMesh, objects, output);
    convertPointFields<tensor>(mesh, pMesh, objects, output);

    if (debug)
    {
        printMemory();
    }
}


void Foam::vtkPVFoam::convertlagrangianFields(vtkMultiBlockDataSet* output)
{
    arrayRange& range = arrayRangelagrangian_;
    const fvMesh& mesh = *meshPtr_;

    wordHashSet selectedFields = getSelected
    (
        reader_->GetlagrangianFieldSelection()
    );

    if (selectedFields.empty())
    {
        return;
    }

    if (debug)
    {
        InfoInFunction << endl;
        printMemory();
    }

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word  cloudName = getPartName(partId);
        const label datasetNo = partDataset_[partId];

        if (!partStatus_[partId] || datasetNo < 0)
        {
            continue;
        }

        // Get the lagrangian fields for this time and this cloud
        // but only keep selected fields
        // the region name is already in the mesh db
        IOobjectList objects
        (
            getObjects
            (
                selectedFields,
                mesh,
                dbPtr_().name(),
                lagrangian::cloud::prefix/cloudName
            )
        );

        if (objects.empty())
        {
            continue;
        }

        if (debug)
        {
            InfoInFunction
                << "converting OpenFOAM lagrangian fields" << nl << "    ";
            forAllConstIter(IOobjectList, objects, iter)
            {
                Info<< "  " << iter()->name()
                    << " == " << iter()->objectPath(false) << nl;
            }
        }

        convertlagrangianFields<label>
        (
            objects, output, datasetNo
        );
        convertlagrangianFields<scalar>
        (
            objects, output, datasetNo
        );
        convertlagrangianFields<vector>
        (
            objects, output, datasetNo
        );
        convertlagrangianFields<sphericalTensor>
        (
            objects, output, datasetNo
        );
        convertlagrangianFields<symmTensor>
        (
            objects, output, datasetNo
        );
        convertlagrangianFields<tensor>
        (
            objects, output, datasetNo
        );
    }

    if (debug)
    {
        printMemory();
    }
}


void Foam::vtkPVFoam::convertLagrangianFields
(
    vtkMultiBlockDataSet* output,
    const PtrList<LagrangianMesh>& LmeshPtrs
)
{
    arrayRange& range = arrayRangeLagrangian_;
    const fvMesh& mesh = *meshPtr_;

    wordHashSet selectedFields = getSelected
    (
        reader_->GetLagrangianFieldSelection()
    );

    if (selectedFields.empty())
    {
        return;
    }

    if (debug)
    {
        InfoInFunction << endl;
        printMemory();
    }

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word  LagrangianName = getPartName(partId);
        const label datasetNo = partDataset_[partId];

        if (!partStatus_[partId] || datasetNo < 0)
        {
            continue;
        }

        // Get the lagrangian fields for this time and this cloud
        // but only keep selected fields
        // the region name is already in the mesh db
        IOobjectList objects
        (
            getObjects
            (
                selectedFields,
                mesh,
                dbPtr_().name(),
                LagrangianMesh::prefix/LagrangianName
            )
        );

        if (objects.empty())
        {
            continue;
        }

        if (debug)
        {
            InfoInFunction
                << "converting OpenFOAM Lagrangian fields" << nl << "    ";
            forAllConstIter(IOobjectList, objects, iter)
            {
                Info<< "  " << iter()->name()
                    << " == " << iter()->objectPath(false) << nl;
            }
        }

        #define CONVERT_LAGRANGIAN_FIELDS(Type, GeoField) \
            convertLagrangianFields<Type, GeoField>       \
            (                                             \
                objects,                                  \
                output,                                   \
                datasetNo,                                \
                LmeshPtrs[datasetNo]                      \
            );
        CONVERT_LAGRANGIAN_FIELDS(label, LagrangianField);
        FOR_ALL_FIELD_TYPES(CONVERT_LAGRANGIAN_FIELDS, LagrangianField);
        CONVERT_LAGRANGIAN_FIELDS(label, LagrangianInternalField);
        FOR_ALL_FIELD_TYPES(CONVERT_LAGRANGIAN_FIELDS, LagrangianInternalField);
        #undef CONVERT_LAGRANGIAN_FIELDS
    }

    if (debug)
    {
        printMemory();
    }
}


// ************************************************************************* //
