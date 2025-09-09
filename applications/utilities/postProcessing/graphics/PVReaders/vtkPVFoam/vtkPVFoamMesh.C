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

// OpenFOAM includes
#include "domainDecomposition.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "fvMeshSubset.H"
#include "lagrangianFieldReconstructor.H"
#include "LagrangianMesh.H"
#include "LagrangianFieldReconstructor.H"
#include "uindirectPrimitivePatch.H"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPVFoam::convertMeshVolume
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    DebugInFunction;

    arrayRange& range = arrayRangeVolume_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0

    const fvMesh& mesh = procMeshesPtr_->completeMesh();

    // Resize for decomposed polyhedra
    regionPolyDecomp_.setSize(range.size());

    // Convert the internalMesh
    // This looks like more than one part, but it isn't
    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word partName = "internalMesh";

        if (!partStatus_[partId]) continue;

        vtkUnstructuredGrid* vtkmesh =
            volumeVTKMesh(mesh, regionPolyDecomp_[datasetNo]);

        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, range, datasetNo, partName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo) ++blockNo;
}


void Foam::vtkPVFoam::convertMeshlagrangian
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    DebugInFunction;

    arrayRange& range = arrayRangelagrangian_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0

    lagrangianReconstructors_.clear();

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word lagrangianName = getPartName(partId);

        if (!partStatus_[partId]) continue;

        autoPtr<lagrangianFieldReconstructor> lreconstructorPtr;

        vtkPolyData* vtkmesh =
            lagrangianVTKMesh(lagrangianName, lreconstructorPtr);
        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, range, datasetNo, lagrangianName);
            vtkmesh->Delete();

            lagrangianReconstructors_.append(lreconstructorPtr.ptr());

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo) ++blockNo;
}


void Foam::vtkPVFoam::convertMeshLagrangian
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    DebugInFunction;

    arrayRange& range = arrayRangeLagrangian_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0

    LagrangianMeshes_.clear();
    LagrangianReconstructors_.clear();

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word LagrangianName = getPartName(partId);

        if (!partStatus_[partId]) continue;

        autoPtr<LagrangianMesh> LmeshPtr;
        autoPtr<LagrangianFieldReconstructor> LreconstructorPtr;

        vtkPolyData* vtkmesh =
            LagrangianVTKMesh(LagrangianName, LmeshPtr, LreconstructorPtr);
        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, range, datasetNo, LagrangianName);
            vtkmesh->Delete();

            LagrangianMeshes_.append(LmeshPtr.ptr());
            LagrangianReconstructors_.append(LreconstructorPtr.ptr());

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo) ++blockNo;
}


void Foam::vtkPVFoam::convertMeshPatches
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangePatches_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0

    const fvMesh& mesh = procMeshesPtr_->completeMesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        if (!partStatus_[partId]) continue;

        const word patchName = getPartName(partId);
        const labelHashSet patchIds =
            patches.patchSet(List<wordRe>(1, wordRe(patchName)));

        DebugInFunction << "Creating VTK mesh for patch(es)[";
        forAllConstIter(labelHashSet, patchIds, iter)
        {
            if (iter != patchIds.begin()) DebugInfo<< ',';
            DebugInfo<< iter.key();
        }
        DebugInfo<< "]: " << patchName << endl;

        vtkPolyData* vtkmesh = nullptr;
        if (patchIds.size() == 1)
        {
            vtkmesh = patchVTKMesh(patches[patchIds.begin().key()]);
        }
        else
        {
            // Patch group. Collect patch faces.
            label sz = 0;
            forAllConstIter(labelHashSet, patchIds, iter)
            {
                sz += patches[iter.key()].size();
            }
            labelList meshFaceLabels(sz);
            sz = 0;
            forAllConstIter(labelHashSet, patchIds, iter)
            {
                const polyPatch& pp = patches[iter.key()];
                forAll(pp, i)
                {
                    meshFaceLabels[sz++] = pp.start()+i;
                }
            }
            UIndirectList<face> fcs(mesh.faces(), meshFaceLabels);
            uindirectPrimitivePatch pp(fcs, mesh.points());

            vtkmesh = patchVTKMesh(pp);
        }

        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, range, datasetNo, patchName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo)
    {
        ++blockNo;
    }
}


void Foam::vtkPVFoam::convertMeshCellZones
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangeCellZones_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0

    const fvMesh& mesh = procMeshesPtr_->completeMesh();

    // Resize for decomposed polyhedra
    zonePolyDecomp_.setSize(range.size());

    if (range.empty()) return;

    const cellZoneList& zMesh = mesh.cellZones();
    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word zoneName = getPartName(partId);
        const label zoneId = zMesh.findIndex(zoneName);

        if (!partStatus_[partId] || zoneId < 0) continue;

        DebugInFunction
            << "Creating VTK mesh for cellZone[" << zoneId << "]: "
            << zoneName << endl;

        fvMeshSubset subsetter(mesh);
        subsetter.setLargeCellSubset(zMesh[zoneId]);

        vtkUnstructuredGrid* vtkmesh = volumeVTKMesh
        (
            subsetter.subMesh(),
            zonePolyDecomp_[datasetNo]
        );

        if (vtkmesh)
        {
            // superCells + addPointCellLabels must contain global cell ids
            inplaceRenumber
            (
                subsetter.cellMap(),
                zonePolyDecomp_[datasetNo].superCells()
            );
            inplaceRenumber
            (
                subsetter.cellMap(),
                zonePolyDecomp_[datasetNo].addPointCellLabels()
            );

            // Copy pointMap as well, otherwise pointFields fail
            zonePolyDecomp_[datasetNo].pointMap() = subsetter.pointMap();

            AddToBlock(output, vtkmesh, range, datasetNo, zoneName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo) ++blockNo;
}


void Foam::vtkPVFoam::convertMeshCellSets
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangeCellSets_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0

    const fvMesh& mesh = procMeshesPtr_->completeMesh();

    // Resize for decomposed polyhedra
    setPolyDecomp_.setSize(range.size());

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word partName = getPartName(partId);

        if (!partStatus_[partId]) continue;

        DebugInFunction
            << "Creating VTK mesh for cellSet: [" << partName << endl;

        const autoPtr<cellSet> cSetPtr =
            reader_->GetDecomposedCase()
          ? procMeshesPtr_->reconstructSet<cellSet>(partName)
          : autoPtr<cellSet>(new cellSet(mesh, partName));

        fvMeshSubset subsetter(mesh);
        subsetter.setLargeCellSubset(cSetPtr());

        vtkUnstructuredGrid* vtkmesh = volumeVTKMesh
        (
            subsetter.subMesh(),
            setPolyDecomp_[datasetNo]
        );

        if (vtkmesh)
        {
            // superCells + addPointCellLabels must contain global cell ids
            inplaceRenumber
            (
                subsetter.cellMap(),
                setPolyDecomp_[datasetNo].superCells()
            );
            inplaceRenumber
            (
                subsetter.cellMap(),
                setPolyDecomp_[datasetNo].addPointCellLabels()
            );

            // Copy pointMap as well, otherwise pointFields fail
            setPolyDecomp_[datasetNo].pointMap() = subsetter.pointMap();

            AddToBlock(output, vtkmesh, range, datasetNo, partName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo) ++blockNo;
}


void Foam::vtkPVFoam::convertMeshFaceZones
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangeFaceZones_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0

    const fvMesh& mesh = procMeshesPtr_->completeMesh();

    if (range.empty()) return;

    const faceZoneList& zMesh = mesh.faceZones();
    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word zoneName = getPartName(partId);
        const label zoneId = zMesh.findIndex(zoneName);

        if (!partStatus_[partId] || zoneId < 0) continue;

        DebugInFunction
            << "Creating VTK mesh for faceZone[" << zoneId << "]: "
            << zoneName << endl;

        vtkPolyData* vtkmesh = patchVTKMesh(zMesh[zoneId].patch());
        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, range, datasetNo, zoneName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo) ++blockNo;
}


void Foam::vtkPVFoam::convertMeshFaceSets
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangeFaceSets_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0

    const fvMesh& mesh = procMeshesPtr_->completeMesh();

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word partName = getPartName(partId);

        if (!partStatus_[partId]) continue;

        DebugInFunction
            << "Creating VTK mesh for faceSet: " << partName << endl;

        const autoPtr<faceSet> fSetPtr =
            reader_->GetDecomposedCase()
          ? procMeshesPtr_->reconstructSet<faceSet>(partName)
          : autoPtr<faceSet>(new faceSet(mesh, partName));

        vtkPolyData* vtkmesh = faceSetVTKMesh(mesh, fSetPtr());
        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, range, datasetNo, partName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo) ++blockNo;
}


void Foam::vtkPVFoam::convertMeshPointZones
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangePointZones_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0

    const fvMesh& mesh = procMeshesPtr_->completeMesh();

    if (range.empty()) return;

    const pointZoneList& zMesh = mesh.pointZones();
    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word zoneName = getPartName(partId);
        const label zoneId = zMesh.findIndex(zoneName);

        if (!partStatus_[partId] || zoneId < 0) continue;

        DebugInFunction
            << "Creating VTK mesh for pointZone[" << zoneId << "]: "
            << zoneName << endl;

        vtkPolyData* vtkmesh = pointZoneVTKMesh(mesh, zMesh[zoneId]);
        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, range, datasetNo, zoneName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo) ++blockNo;
}



void Foam::vtkPVFoam::convertMeshPointSets
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangePointSets_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0

    const fvMesh& mesh = procMeshesPtr_->completeMesh();

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        word partName = getPartName(partId);

        if (!partStatus_[partId]) continue;

        DebugInFunction
            << "Creating VTK mesh for pointSet: " << partName << endl;

        const autoPtr<pointSet> pSetPtr =
            reader_->GetDecomposedCase()
          ? procMeshesPtr_->reconstructSet<pointSet>(partName)
          : autoPtr<pointSet>(new pointSet(mesh, partName));

        vtkPolyData* vtkmesh = pointSetVTKMesh(mesh, pSetPtr());
        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, range, datasetNo, partName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo) ++blockNo;
}


// ************************************************************************* //
