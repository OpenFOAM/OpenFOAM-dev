/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "blockMesh.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Pair<Foam::scalar> Foam::blockMesh::xCellSizes
(
    const block& b,
    const pointField& blockPoints,
    const label j,
    const label k
) const
{
    return Pair<scalar>
    (
        mag
        (
            blockPoints[b.pointLabel(1, j, k)]
          - blockPoints[b.pointLabel(0, j, k)]
        ),

        mag
        (
            blockPoints[b.pointLabel(b.density().x() - 1, j, k)]
          - blockPoints[b.pointLabel(b.density().x(), j, k)]
        )
    );
}


Foam::Pair<Foam::scalar> Foam::blockMesh::yCellSizes
(
    const block& b,
    const pointField& blockPoints,
    const label i,
    const label k
) const
{
    return Pair<scalar>
    (
        mag
        (
            blockPoints[b.pointLabel(i, 0, k)]
          - blockPoints[b.pointLabel(i, 1, k)]
        ),

        mag
        (
            blockPoints[b.pointLabel(i, b.density().y() - 1, k)]
          - blockPoints[b.pointLabel(i, b.density().y(), k)]
        )
    );
}


Foam::Pair<Foam::scalar> Foam::blockMesh::zCellSizes
(
    const block& b,
    const pointField& blockPoints,
    const label i,
    const label j
) const
{
    return Pair<scalar>
    (
        mag
        (
            blockPoints[b.pointLabel(i, j, 0)]
          - blockPoints[b.pointLabel(i, j, 1)]
        ),

        mag
        (
            blockPoints[b.pointLabel(i, j, b.density().z() - 1)]
          - blockPoints[b.pointLabel(i, j, b.density().z())]
        )
    );
}


void Foam::blockMesh::printCellSizeRange
(
    const Pair<scalar>& cellSizes
) const
{
    if (cellSizes.first() != cellSizes.second())
    {
        Info<< scaleFactor_*cellSizes.first() << " .. "
            << scaleFactor_*cellSizes.second();
    }
    else
    {
        Info<< scaleFactor_*cellSizes.first();
    }
}


void Foam::blockMesh::printCellSizeRanges
(
    const direction d,
    const FixedList<Pair<scalar>, 4>& cellSizes
) const
{
    static const char dNames[3] = {'i', 'j', 'k'};

    const scalar d0 = max
    (
        max(cellSizes[0].first(), cellSizes[0].second()),
        small
    );

    bool uniform = true;
    forAll(cellSizes, i)
    {
        uniform = uniform
           && mag(cellSizes[i].first() - cellSizes[0].first())/d0 < 1e-5
           && mag(cellSizes[i].second() - cellSizes[0].second())/d0 < 1e-5;
    }

    Info<< "        " << dNames[d] << " : ";
    if (uniform)
    {
        printCellSizeRange(cellSizes[0]);
    }
    else
    {
        forAll(cellSizes, i)
        {
            printCellSizeRange(cellSizes[i]);
            Info<< " ";
        }
    }
    Info<< endl;
}


void Foam::blockMesh::createPoints() const
{
    const blockList& blocks = *this;

    if (verboseOutput)
    {
        Info<< "Creating points with scale " << scaleFactor_ << endl;
    }

    points_.setSize(nPoints_);

    forAll(blocks, blocki)
    {
        const pointField& blockPoints = blocks[blocki].points();

        if (verboseOutput)
        {
            Info<< "    Block " << blocki << " cell size :" << nl;

            printCellSizeRanges
            (
                0,
                {
                    xCellSizes(blocks[blocki], blockPoints, 0, 0),
                    xCellSizes(blocks[blocki], blockPoints, 1, 0),
                    xCellSizes(blocks[blocki], blockPoints, 0, 1),
                    xCellSizes(blocks[blocki], blockPoints, 1, 1)
                }
            );

            printCellSizeRanges
            (
                1,
                {
                    yCellSizes(blocks[blocki], blockPoints, 0, 0),
                    yCellSizes(blocks[blocki], blockPoints, 1, 0),
                    yCellSizes(blocks[blocki], blockPoints, 0, 1),
                    yCellSizes(blocks[blocki], blockPoints, 1, 1)
                }
            );

            printCellSizeRanges
            (
                2,
                {
                    zCellSizes(blocks[blocki], blockPoints, 0, 0),
                    zCellSizes(blocks[blocki], blockPoints, 1, 0),
                    zCellSizes(blocks[blocki], blockPoints, 0, 1),
                    zCellSizes(blocks[blocki], blockPoints, 1, 1)
                }
            );
        }

        forAll(blockPoints, blockPointi)
        {
            points_
            [
                mergeList_
                [
                    blockOffsets_[blocki] + blockPointi
                ]
            ] = scaleFactor_ * blockPoints[blockPointi];
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

    cells_.setSize(nCells_);

    label cellLabel = 0;

    forAll(blocks, blocki)
    {
        const List<FixedList<label, 8>> blockCells(blocks[blocki].cells());

        forAll(blockCells, blockCelli)
        {
            labelList cellPoints(blockCells[blockCelli].size());

            forAll(cellPoints, cellPointi)
            {
                cellPoints[cellPointi] =
                    mergeList_
                    [
                        blockCells[blockCelli][cellPointi]
                      + blockOffsets_[blocki]
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
        const label blocki = blockLabels[patchTopologyFaceLabel];

        faceList blockFaces = blocks[blocki].blockShape().faces();

        forAll(blockFaces, blockFaceLabel)
        {
            if
            (
                blockFaces[blockFaceLabel]
             == patchTopologyFaces[patchTopologyFaceLabel]
            )
            {
                nFaces +=
                    blocks[blocki].boundaryPatches()[blockFaceLabel].size();
            }
        }
    }


    faceList patchFaces(nFaces);
    face quadFace(4);
    label faceLabel = 0;

    forAll(patchTopologyFaces, patchTopologyFaceLabel)
    {
        const label blocki = blockLabels[patchTopologyFaceLabel];

        faceList blockFaces = blocks[blocki].blockShape().faces();

        forAll(blockFaces, blockFaceLabel)
        {
            if
            (
                blockFaces[blockFaceLabel]
             == patchTopologyFaces[patchTopologyFaceLabel]
            )
            {
                const List<FixedList<label, 4>>& blockPatchFaces =
                    blocks[blocki].boundaryPatches()[blockFaceLabel];

                forAll(blockPatchFaces, blockFaceLabel)
                {
                    // Lookup the face points
                    // and collapse duplicate point labels

                    quadFace[0] =
                        mergeList_
                        [
                            blockPatchFaces[blockFaceLabel][0]
                          + blockOffsets_[blocki]
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
                              + blockOffsets_[blocki]
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

    patches_.setSize(topoPatches.size());

    forAll(topoPatches, patchi)
    {
        patches_[patchi] = createPatchFaces(topoPatches[patchi]);
    }
}


// ************************************************************************* //
