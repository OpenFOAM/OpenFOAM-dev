/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "error.H"
#include "blockMesh.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::blockMesh::createPoints() const
{
    const blockList& blocks = *this;

    if (verboseOutput)
    {
        Info<< "Creating points with scale " << scaleFactor_ << endl;
    }

    //
    // generate points
    //
    points_.clear();
    points_.setSize(nPoints_);

    forAll(blocks, blockI)
    {
        const pointField& blockPoints = blocks[blockI].points();

        if (verboseOutput)
        {
            const Vector<label>& density = blocks[blockI].meshDensity();

            label v0 = blocks[blockI].vtxLabel(0, 0, 0);
            label vi1 = blocks[blockI].vtxLabel(1, 0, 0);
            scalar diStart = mag(blockPoints[vi1]-blockPoints[v0]);

            label vinM1 = blocks[blockI].vtxLabel(density.x()-1, 0, 0);
            label vin = blocks[blockI].vtxLabel(density.x(), 0, 0);
            scalar diFinal = mag(blockPoints[vin]-blockPoints[vinM1]);

            label vj1 = blocks[blockI].vtxLabel(0, 1, 0);
            scalar djStart = mag(blockPoints[vj1]-blockPoints[v0]);
            label vjnM1 = blocks[blockI].vtxLabel(0, density.y()-1, 0);
            label vjn = blocks[blockI].vtxLabel(0, density.y(), 0);
            scalar djFinal = mag(blockPoints[vjn]-blockPoints[vjnM1]);

            label vk1 = blocks[blockI].vtxLabel(0, 0, 1);
            scalar dkStart = mag(blockPoints[vk1]-blockPoints[v0]);
            label vknM1 = blocks[blockI].vtxLabel(0, 0, density.z()-1);
            label vkn = blocks[blockI].vtxLabel(0, 0, density.z());
            scalar dkFinal = mag(blockPoints[vkn]-blockPoints[vknM1]);

            Info<< "    Block " << blockI << " cell size :" << nl
                << "        i : " << scaleFactor_*diStart << " .. "
                << scaleFactor_*diFinal << nl
                << "        j : " << scaleFactor_*djStart << " .. "
                << scaleFactor_*djFinal << nl
                << "        k : " << scaleFactor_*dkStart << " .. "
                << scaleFactor_*dkFinal << nl
                << endl;
        }

        forAll(blockPoints, blockPointI)
        {
            points_
            [
                mergeList_
                [
                    blockOffsets_[blockI] + blockPointI
                ]
            ] = scaleFactor_ * blockPoints[blockPointI];
        }
    }
}


void Foam::blockMesh::createCells() const
{
    const blockList& blocks = *this;
    const cellModel& hex = *(cellModeller::lookup("hex"));

    if (verboseOutput)
    {
        Info<< "Creating cells" << endl;
    }

    //
    // generate cells
    //
    cells_.clear();
    cells_.setSize(nCells_);

    label cellLabel = 0;

    forAll(blocks, blockI)
    {
        const labelListList& blockCells = blocks[blockI].cells();

        forAll(blockCells, blockCellI)
        {
            labelList cellPoints(blockCells[blockCellI].size());

            forAll(cellPoints, cellPointI)
            {
                cellPoints[cellPointI] =
                    mergeList_
                    [
                        blockCells[blockCellI][cellPointI]
                      + blockOffsets_[blockI]
                    ];
            }

            // Construct collapsed cell and add to list
            cells_[cellLabel] = cellShape(hex, cellPoints, true);

            cellLabel++;
        }
    }
}


Foam::faceList Foam::blockMesh::createPatchFaces
(
    const polyPatch& patchTopologyFaces
) const
{
    const blockList& blocks = *this;

    labelList blockLabels = patchTopologyFaces.polyPatch::faceCells();

    label nFaces = 0;

    forAll(patchTopologyFaces, patchTopologyFaceLabel)
    {
        const label blockI = blockLabels[patchTopologyFaceLabel];

        faceList blockFaces = blocks[blockI].blockShape().faces();

        forAll(blockFaces, blockFaceLabel)
        {
            if
            (
                blockFaces[blockFaceLabel]
             == patchTopologyFaces[patchTopologyFaceLabel]
            )
            {
                nFaces +=
                    blocks[blockI].boundaryPatches()[blockFaceLabel].size();
            }
        }
    }


    faceList patchFaces(nFaces);
    face quadFace(4);
    label faceLabel = 0;

    forAll(patchTopologyFaces, patchTopologyFaceLabel)
    {
        const label blockI = blockLabels[patchTopologyFaceLabel];

        faceList blockFaces = blocks[blockI].blockShape().faces();

        forAll(blockFaces, blockFaceLabel)
        {
            if
            (
                blockFaces[blockFaceLabel]
                == patchTopologyFaces[patchTopologyFaceLabel]
            )
            {
                const labelListList& blockPatchFaces =
                    blocks[blockI].boundaryPatches()[blockFaceLabel];

                forAll(blockPatchFaces, blockFaceLabel)
                {
                    // Lookup the face points
                    // and collapse duplicate point labels

                    quadFace[0] =
                        mergeList_
                        [
                            blockPatchFaces[blockFaceLabel][0]
                          + blockOffsets_[blockI]
                        ];

                    label nUnique = 1;

                    for
                    (
                        label facePointLabel = 1;
                        facePointLabel < 4;
                        facePointLabel++
                    )
                    {
                        quadFace[nUnique] =
                            mergeList_
                            [
                                blockPatchFaces[blockFaceLabel][facePointLabel]
                              + blockOffsets_[blockI]
                            ];

                        if (quadFace[nUnique] != quadFace[nUnique-1])
                        {
                            nUnique++;
                        }
                    }

                    if (quadFace[nUnique-1] == quadFace[0])
                    {
                        nUnique--;
                    }

                    if (nUnique == 4)
                    {
                        patchFaces[faceLabel++] = quadFace;
                    }
                    else if (nUnique == 3)
                    {
                        patchFaces[faceLabel++] = face
                        (
                            labelList::subList(quadFace, 3)
                        );
                    }
                    // else the face has collapsed to an edge or point
                }
            }
        }
    }

    patchFaces.setSize(faceLabel);

    return patchFaces;
}


void Foam::blockMesh::createPatches() const
{
    const polyPatchList& topoPatches = topology().boundaryMesh();

    if (verboseOutput)
    {
        Info<< "Creating patches" << endl;
    }

    //
    // generate points
    //

    patches_.clear();
    patches_.setSize(topoPatches.size());

    forAll(topoPatches, patchI)
    {
        patches_[patchI] = createPatchFaces(topoPatches[patchI]);
    }

}


void Foam::blockMesh::clearGeom()
{
    blockList& blocks = *this;

    forAll(blocks, blockI)
    {
        blocks[blockI].clearGeom();
    }
}

// ************************************************************************* //
