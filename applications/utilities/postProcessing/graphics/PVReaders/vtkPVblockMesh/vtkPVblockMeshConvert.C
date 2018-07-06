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

\*---------------------------------------------------------------------------*/

#include "vtkPVblockMesh.H"
#include "vtkPVblockMeshReader.h"

// OpenFOAM includes
#include "blockMesh.H"
#include "Time.H"

#include "vtkOpenFOAMPoints.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPVblockMesh::convertMeshBlocks
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    vtkDataArraySelection* selection = reader_->GetBlockSelection();
    arrayRange& range = arrayRangeBlocks_;
    range.block(blockNo);   // set output block
    label datasetNo = 0;       // restart at dataset 0

    const blockMesh& blkMesh = *meshPtr_;
    const Foam::pointField& blockPoints = blkMesh.vertices();

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVblockMesh::convertMeshBlocks" << endl;
    }

    int blockI = 0;
    const scalar scaleFactor = blkMesh.scaleFactor();

    for
    (
        int partId = range.start();
        partId < range.end();
        ++partId, ++blockI
    )
    {
        if (!blockStatus_[partId])
        {
            continue;
        }

        const blockDescriptor& blockDef = blkMesh[blockI];

        vtkUnstructuredGrid* vtkmesh = vtkUnstructuredGrid::New();

        // Convert OpenFOAM mesh vertices to VTK
        vtkPoints *vtkpoints = vtkPoints::New();
        vtkpoints->Allocate( blockDef.nPoints() );
        const labelList& blockLabels = blockDef.blockShape();

        vtkmesh->Allocate(1);
        vtkIdType nodeIds[8];

        forAll(blockLabels, ptI)
        {
            vtkInsertNextOpenFOAMPoint
            (
                vtkpoints,
                blockPoints[blockLabels[ptI]],
                scaleFactor
            );

            nodeIds[ptI] = ptI;
        }

        vtkmesh->InsertNextCell
        (
            VTK_HEXAHEDRON,
            8,
            nodeIds
        );

        vtkmesh->SetPoints(vtkpoints);
        vtkpoints->Delete();

        AddToBlock
        (
            output, vtkmesh, range, datasetNo,
            selection->GetArrayName(partId)
        );

        vtkmesh->Delete();
        datasetNo++;
    }


    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPVblockMesh::convertMeshBlocks" << endl;
    }
}


void Foam::vtkPVblockMesh::convertMeshEdges
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    vtkDataArraySelection* selection = reader_->GetCurvedEdgesSelection();
    arrayRange& range = arrayRangeEdges_;

    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0

    const blockMesh& blkMesh = *meshPtr_;
    const blockEdgeList& edges = blkMesh.edges();

    int edgeI = 0;
    const scalar scaleFactor = blkMesh.scaleFactor();

    for
    (
        int partId = range.start();
        partId < range.end();
        ++partId, ++edgeI
    )
    {
        if (!edgeStatus_[partId])
        {
            continue;
        }

        // search each block
        forAll(blkMesh, blockI)
        {
            const blockDescriptor& blockDef = blkMesh[blockI];

            edgeList blkEdges = blockDef.blockShape().edges();

            // List of edge point and weighting factors
            pointField edgesPoints[12];
            scalarList edgesWeights[12];
            blockDef.edgesPointsWeights(edgesPoints, edgesWeights);

            // find the corresponding edge within the block
            label foundEdgeI = -1;
            forAll(blkEdges, blkEdgeI)
            {
                if (edges[edgeI].compare(blkEdges[blkEdgeI]))
                {
                    foundEdgeI = blkEdgeI;
                    break;
                }
            }

            if (foundEdgeI != -1)
            {
                const List<point>& edgePoints = edgesPoints[foundEdgeI];

                vtkPolyData* vtkmesh = vtkPolyData::New();
                vtkPoints* vtkpoints = vtkPoints::New();

                vtkpoints->Allocate( edgePoints.size() );
                vtkmesh->Allocate(1);

                vtkIdType pointIds[edgePoints.size()];
                forAll(edgePoints, ptI)
                {
                    vtkInsertNextOpenFOAMPoint
                    (
                        vtkpoints,
                        edgePoints[ptI],
                        scaleFactor
                    );
                    pointIds[ptI] = ptI;
                }

                vtkmesh->InsertNextCell
                (
                    VTK_POLY_LINE,
                    edgePoints.size(),
                    pointIds
                );

                vtkmesh->SetPoints(vtkpoints);
                vtkpoints->Delete();

                AddToBlock
                (
                    output, vtkmesh, range, datasetNo,
                    selection->GetArrayName(partId)
                );

                vtkmesh->Delete();
                datasetNo++;

                break;
            }
        }
    }


    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPVblockMesh::convertMeshEdges" << endl;
    }

}


void Foam::vtkPVblockMesh::convertMeshCorners
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangeCorners_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0

    const pointField& blockPoints = meshPtr_->vertices();
    const scalar& scaleFactor = meshPtr_->scaleFactor();

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVblockMesh::convertMeshCorners" << endl;
    }

    if (true)  // or some flag or other condition
    {
        vtkPolyData* vtkmesh = vtkPolyData::New();
        vtkPoints* vtkpoints = vtkPoints::New();
        vtkCellArray* vtkcells = vtkCellArray::New();

        vtkpoints->Allocate( blockPoints.size() );
        vtkcells->Allocate( blockPoints.size() );

        vtkIdType pointId = 0;
        forAll(blockPoints, ptI)
        {
            vtkInsertNextOpenFOAMPoint
            (
                vtkpoints,
                blockPoints[ptI],
                scaleFactor
            );

            vtkcells->InsertNextCell(1, &pointId);
            pointId++;
        }

        vtkmesh->SetPoints(vtkpoints);
        vtkpoints->Delete();

        vtkmesh->SetVerts(vtkcells);
        vtkcells->Delete();

        AddToBlock
        (
            output, vtkmesh, range, datasetNo,
            arrayRangeCorners_.name()
        );
        vtkmesh->Delete();

        datasetNo++;
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPVblockMesh::convertMeshCorners" << endl;
    }
}


// ************************************************************************* //
