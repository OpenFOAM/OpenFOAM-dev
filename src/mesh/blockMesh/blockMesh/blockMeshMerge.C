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

#include "blockMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::blockMesh::calcMergeInfo()
{
    const blockList& blocks = *this;

    if (verboseOutput)
    {
        Info<< "Creating block offsets" << endl;
    }

    blockOffsets_.setSize(blocks.size());

    nPoints_ = 0;
    nCells_  = 0;

    forAll(blocks, blocki)
    {
        blockOffsets_[blocki] = nPoints_;

        nPoints_ += blocks[blocki].nPoints();
        nCells_  += blocks[blocki].nCells();
    }


    if (verboseOutput)
    {
        Info<< "Creating merge list " << flush;
    }

    // set unused to -1
    mergeList_.setSize(nPoints_);
    mergeList_ = -1;


    const pointField& blockPoints = topology().points();
    const cellList& blockCells = topology().cells();
    const faceList& blockFaces = topology().faces();
    const labelList& faceOwnerBlocks = topology().faceOwner();

    // For efficiency, create merge pairs in the first pass
    labelListListList glueMergePairs(blockFaces.size());

    const labelList& faceNeighbourBlocks = topology().faceNeighbour();


    forAll(blockFaces, blockFaceLabel)
    {
        label blockPlabel = faceOwnerBlocks[blockFaceLabel];
        const pointField& blockPpoints = blocks[blockPlabel].points();
        const labelList& blockPfaces = blockCells[blockPlabel];

        bool foundFace = false;
        label blockPfaceLabel;
        for
        (
            blockPfaceLabel = 0;
            blockPfaceLabel < blockPfaces.size();
            blockPfaceLabel++
        )
        {
            if (blockPfaces[blockPfaceLabel] == blockFaceLabel)
            {
                foundFace = true;
                break;
            }
        }

        if (!foundFace)
        {
            FatalErrorInFunction
                << "Cannot find merge face for block " << blockPlabel
                << exit(FatalError);
        }

        const List<FixedList<label, 4>>& blockPfaceFaces =
            blocks[blockPlabel].boundaryPatches()[blockPfaceLabel];

        labelListList& curPairs = glueMergePairs[blockFaceLabel];
        curPairs.setSize(blockPfaceFaces.size());

        // Calculate sqr of the merge tolerance as 1/10th of the min sqr
        // point to point distance on the block face.
        // At the same time merge collated points on the block's faces
        // (removes boundary poles etc.)
        // Collated points detected by initially taking a constant factor of
        // the size of the block.

        boundBox bb(blockCells[blockPlabel].points(blockFaces, blockPoints));
        const scalar mergeSqrDist = magSqr(10*small*bb.span());

        // This is an N^2 algorithm

        scalar sqrMergeTol = great;

        forAll(blockPfaceFaces, blockPfaceFaceLabel)
        {
            const FixedList<label, 4>& blockPfaceFacePoints
                = blockPfaceFaces[blockPfaceFaceLabel];

            forAll(blockPfaceFacePoints, blockPfaceFacePointLabel)
            {
                forAll(blockPfaceFacePoints, blockPfaceFacePointLabel2)
                {
                    if (blockPfaceFacePointLabel != blockPfaceFacePointLabel2)
                    {
                        scalar magSqrDist = magSqr
                        (
                            blockPpoints[blockPfaceFacePoints
                            [blockPfaceFacePointLabel]]
                            - blockPpoints[blockPfaceFacePoints
                            [blockPfaceFacePointLabel2]]
                        );

                        if (magSqrDist < mergeSqrDist)
                        {
                            label PpointLabel =
                                blockPfaceFacePoints[blockPfaceFacePointLabel]
                              + blockOffsets_[blockPlabel];

                            label PpointLabel2 =
                                blockPfaceFacePoints[blockPfaceFacePointLabel2]
                              + blockOffsets_[blockPlabel];

                            label minPP2 = min(PpointLabel, PpointLabel2);

                            if (mergeList_[PpointLabel] != -1)
                            {
                                minPP2 = min(minPP2, mergeList_[PpointLabel]);
                            }

                            if (mergeList_[PpointLabel2] != -1)
                            {
                                minPP2 = min(minPP2, mergeList_[PpointLabel2]);
                            }

                            mergeList_[PpointLabel] = mergeList_[PpointLabel2]
                                = minPP2;
                        }
                        else
                        {
                            sqrMergeTol = min(sqrMergeTol, magSqrDist);
                        }
                    }
                }
            }
        }

        sqrMergeTol /= 10.0;


        if (topology().isInternalFace(blockFaceLabel))
        {
        label blockNlabel = faceNeighbourBlocks[blockFaceLabel];
        const pointField& blockNpoints = blocks[blockNlabel].points();
        const labelList& blockNfaces = blockCells[blockNlabel];

        foundFace = false;
        label blockNfaceLabel;
        for
        (
            blockNfaceLabel = 0;
            blockNfaceLabel < blockNfaces.size();
            blockNfaceLabel++
        )
        {
            if
            (
                blockFaces[blockNfaces[blockNfaceLabel]]
             == blockFaces[blockFaceLabel]
            )
            {
                foundFace = true;
                break;
            }
        }

        if (!foundFace)
        {
            FatalErrorInFunction
                << "Cannot find merge face for block " << blockNlabel
                << exit(FatalError);
        }

        const List<FixedList<label, 4>>& blockNfaceFaces =
            blocks[blockNlabel].boundaryPatches()[blockNfaceLabel];

        if
        (
            checkFaceCorrespondence_
         && blockPfaceFaces.size() != blockNfaceFaces.size()
        )
        {
            FatalErrorInFunction
                << "Inconsistent number of faces between block pair "
                << blockPlabel << " and " << blockNlabel
                << exit(FatalError);
        }


        // N-squared point search over all points of all faces of
        // master block over all point of all faces of slave block
        forAll(blockPfaceFaces, blockPfaceFaceLabel)
        {
            const FixedList<label, 4>& blockPfaceFacePoints
                = blockPfaceFaces[blockPfaceFaceLabel];

            labelList& cp = curPairs[blockPfaceFaceLabel];
            cp.setSize(blockPfaceFacePoints.size());
            cp = -1;

            forAll(blockPfaceFacePoints, blockPfaceFacePointLabel)
            {
                forAll(blockNfaceFaces, blockNfaceFaceLabel)
                {
                    const FixedList<label, 4>& blockNfaceFacePoints
                        = blockNfaceFaces[blockNfaceFaceLabel];

                    forAll(blockNfaceFacePoints, blockNfaceFacePointLabel)
                    {
                        if
                        (
                            magSqr
                            (
                                blockPpoints
                                [blockPfaceFacePoints[blockPfaceFacePointLabel]]
                              - blockNpoints
                                [blockNfaceFacePoints[blockNfaceFacePointLabel]]
                            ) < sqrMergeTol
                        )
                        {
                            // Found a new pair

                            cp[blockPfaceFacePointLabel] =
                                blockNfaceFacePoints[blockNfaceFacePointLabel];

                            label PpointLabel =
                                blockPfaceFacePoints[blockPfaceFacePointLabel]
                              + blockOffsets_[blockPlabel];

                            label NpointLabel =
                                blockNfaceFacePoints[blockNfaceFacePointLabel]
                              + blockOffsets_[blockNlabel];

                            label minPN = min(PpointLabel, NpointLabel);

                            if (mergeList_[PpointLabel] != -1)
                            {
                                minPN = min(minPN, mergeList_[PpointLabel]);
                            }

                            if (mergeList_[NpointLabel] != -1)
                            {
                                minPN = min(minPN, mergeList_[NpointLabel]);
                            }

                            mergeList_[PpointLabel] = mergeList_[NpointLabel]
                                = minPN;
                        }
                    }
                }
            }
            forAll(blockPfaceFacePoints, blockPfaceFacePointLabel)
            {
                if (cp[blockPfaceFacePointLabel] == -1)
                {
                    FatalErrorInFunction
                        << "Inconsistent point locations between block pair "
                        << blockPlabel << " and " << blockNlabel << nl
                        << "    probably due to inconsistent grading."
                        << exit(FatalError);
                }
            }
        }
        }
    }


    const faceList::subList blockinternalFaces
    (
        blockFaces,
        topology().nInternalFaces()
    );

    bool changedPointMerge = false;
    label nPasses = 0;

    do
    {
        changedPointMerge = false;
        nPasses++;

        forAll(blockinternalFaces, blockFaceLabel)
        {
            label blockPlabel = faceOwnerBlocks[blockFaceLabel];
            label blockNlabel = faceNeighbourBlocks[blockFaceLabel];

            const labelList& blockPfaces = blockCells[blockPlabel];
            const labelList& blockNfaces = blockCells[blockNlabel];

            const labelListList& curPairs = glueMergePairs[blockFaceLabel];

            label blockPfaceLabel;
            for
            (
                blockPfaceLabel = 0;
                blockPfaceLabel < blockPfaces.size();
                blockPfaceLabel++
            )
            {
                if
                (
                    blockFaces[blockPfaces[blockPfaceLabel]]
                 == blockinternalFaces[blockFaceLabel]
                )
                {
                    break;
                }
            }


            label blockNfaceLabel;
            for
            (
                blockNfaceLabel = 0;
                blockNfaceLabel < blockNfaces.size();
                blockNfaceLabel++
            )
            {
                if
                (
                    blockFaces[blockNfaces[blockNfaceLabel]]
                 == blockinternalFaces[blockFaceLabel]
                )
                {
                    break;
                }
            }


            const List<FixedList<label, 4>>& blockPfaceFaces =
                blocks[blockPlabel].boundaryPatches()[blockPfaceLabel];

            forAll(blockPfaceFaces, blockPfaceFaceLabel)
            {
                const FixedList<label, 4>& blockPfaceFacePoints
                    = blockPfaceFaces[blockPfaceFaceLabel];

                const labelList& cp = curPairs[blockPfaceFaceLabel];

                forAll(blockPfaceFacePoints, blockPfaceFacePointLabel)
                {
                    label PpointLabel =
                        blockPfaceFacePoints[blockPfaceFacePointLabel]
                      + blockOffsets_[blockPlabel];

                    label NpointLabel =
                        cp[blockPfaceFacePointLabel]
                      + blockOffsets_[blockNlabel];

                    if
                    (
                        mergeList_[PpointLabel]
                     != mergeList_[NpointLabel]
                    )
                    {
                        changedPointMerge = true;

                        mergeList_[PpointLabel]
                      = mergeList_[NpointLabel]
                      = min
                        (
                            mergeList_[PpointLabel],
                            mergeList_[NpointLabel]
                        );
                    }
                }
            }
        }
        if (verboseOutput)
        {
            Info<< "." << flush;
        }

        if (nPasses > 100)
        {
            FatalErrorInFunction
                << "Point merging failed after max number of passes."
                << exit(FatalError);
        }
    }
    while (changedPointMerge);

    if (verboseOutput)
    {
        Info<< endl;
    }

    forAll(blockinternalFaces, blockFaceLabel)
    {
        label blockPlabel = faceOwnerBlocks[blockFaceLabel];
        label blockNlabel = faceNeighbourBlocks[blockFaceLabel];

        const labelList& blockPfaces = blockCells[blockPlabel];
        const labelList& blockNfaces = blockCells[blockNlabel];

        const pointField& blockPpoints = blocks[blockPlabel].points();
        const pointField& blockNpoints = blocks[blockNlabel].points();

        bool foundFace = false;
        label blockPfaceLabel;
        for
        (
            blockPfaceLabel = 0;
            blockPfaceLabel < blockPfaces.size();
            blockPfaceLabel++
        )
        {
            if
            (
                blockFaces[blockPfaces[blockPfaceLabel]]
             == blockinternalFaces[blockFaceLabel]
            )
            {
                foundFace = true;
                break;
            }
        }

        if (!foundFace)
        {
            FatalErrorInFunction
                << "Cannot find merge face for block " << blockPlabel
                << exit(FatalError);
        }

        foundFace = false;
        label blockNfaceLabel;
        for
        (
            blockNfaceLabel = 0;
            blockNfaceLabel < blockNfaces.size();
            blockNfaceLabel++
        )
        {
            if
            (
                blockFaces[blockNfaces[blockNfaceLabel]]
             == blockinternalFaces[blockFaceLabel]
            )
            {
                foundFace = true;
                break;
            }
        }

        if (!foundFace)
        {
            FatalErrorInFunction
                << "Cannot find merge face for block " << blockNlabel
                << exit(FatalError);
        }

        const List<FixedList<label, 4>>& blockPfaceFaces =
            blocks[blockPlabel].boundaryPatches()[blockPfaceLabel];

        const List<FixedList<label, 4>>& blockNfaceFaces =
            blocks[blockNlabel].boundaryPatches()[blockNfaceLabel];

        forAll(blockPfaceFaces, blockPfaceFaceLabel)
        {
            const FixedList<label, 4>& blockPfaceFacePoints
                = blockPfaceFaces[blockPfaceFaceLabel];

            forAll(blockPfaceFacePoints, blockPfaceFacePointLabel)
            {
                label PpointLabel =
                    blockPfaceFacePoints[blockPfaceFacePointLabel]
                  + blockOffsets_[blockPlabel];

                if (mergeList_[PpointLabel] == -1)
                {
                    FatalErrorInFunction
                        << "Unable to merge point "
                        << blockPfaceFacePointLabel
                        << ' ' << blockPpoints[blockPfaceFacePointLabel]
                        << " of face "
                        << blockPfaceLabel
                        << " of block "
                        << blockPlabel
                        << exit(FatalError);
                }
            }
        }

        forAll(blockNfaceFaces, blockNfaceFaceLabel)
        {
            const FixedList<label, 4>& blockNfaceFacePoints
                = blockNfaceFaces[blockNfaceFaceLabel];

            forAll(blockNfaceFacePoints, blockNfaceFacePointLabel)
            {
                label NpointLabel =
                    blockNfaceFacePoints[blockNfaceFacePointLabel]
                  + blockOffsets_[blockNlabel];

                if (mergeList_[NpointLabel] == -1)
                {
                    FatalErrorInFunction
                        << "unable to merge point "
                        << blockNfaceFacePointLabel
                        << ' ' << blockNpoints[blockNfaceFacePointLabel]
                        << " of face "
                        << blockNfaceLabel
                        << " of block "
                        << blockNlabel
                        << exit(FatalError);
                }
            }
        }
    }


    // Sort merge list to return new point label (in new shorter list)
    // given old point label
    label newPointLabel = 0;

    forAll(mergeList_, pointLabel)
    {
        if (mergeList_[pointLabel] > pointLabel)
        {
            FatalErrorInFunction
                << "Merge list contains point index out of range"
                << exit(FatalError);
        }

        if
        (
            mergeList_[pointLabel] == -1
         || mergeList_[pointLabel] == pointLabel
        )
        {
            mergeList_[pointLabel] = newPointLabel;
            newPointLabel++;
        }
        else
        {
            mergeList_[pointLabel] = mergeList_[mergeList_[pointLabel]];
        }
    }

    nPoints_ = newPointLabel;
}


// ************************************************************************* //
