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

Application
    cfx4ToFoam

Description
    Converts a CFX 4 mesh to OpenFOAM format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IFstream.H"
#include "hexBlock.H"
#include "polyMesh.H"
#include "wallPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "preservePatchTypes.H"
#include "cellShape.H"
#include "cellModeller.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("CFX geom file");
    argList::addOption
    (
        "scale",
        "factor",
        "geometry scaling factor - default is 1"
    );

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }

    const scalar scaleFactor = args.optionLookupOrDefault("scale", 1.0);

#   include "createTime.H"

    IFstream cfxFile(args[1]);

    // Read the cfx information using a fixed format reader.
    // Comments in the file are in C++ style, so the stream parser will remove
    // them with no intervention
    label nblock, npatch, nglue, nelem, npoint;

    cfxFile >> nblock >> npatch >> nglue >> nelem >> npoint;

    Info<< "Reading blocks" << endl;

    PtrList<hexBlock> blocks(nblock);

    {
        word blockName;
        label nx, ny, nz;

        forAll(blocks, blockI)
        {
            cfxFile >> blockName;
            cfxFile >> nx >> ny >> nz;

            blocks.set(blockI, new hexBlock(nx, ny, nz));
        }
    }

    Info<< "Reading patch definitions" << endl;

    wordList cfxPatchTypes(npatch);
    wordList cfxPatchNames(npatch);
    labelList patchMasterBlocks(npatch);
    labelList patchDirections(npatch);
    labelListList patchRanges(npatch);

    {
        label no, blkNo, patchLabel;

        forAll(cfxPatchTypes, patchI)
        {
            // Grab patch type and name
            cfxFile >> cfxPatchTypes[patchI] >> cfxPatchNames[patchI] >> no;

            // Grab patch range
            patchRanges[patchI].setSize(6);
            labelList& curRange = patchRanges[patchI];

            forAll(curRange, rI)
            {
                cfxFile >> curRange[rI];
            }

            // Grab patch direction and master block ID
            // Note: direc is the direction, from the cfx manual
            // 0 = solid (3-D patch),
            // 1 = high i, 2 = high j, 3 = high k
            // 4 = low i, 5 = low j, 6 = low k
            cfxFile >> patchDirections[patchI] >> blkNo >> patchLabel;

            patchMasterBlocks[patchI] = blkNo - 1;
        }
    }

    Info<< "Reading block glueing information" << endl;

    labelList glueMasterPatches(nglue, -1);
    labelList glueSlavePatches(nglue, -1);

    {
        label masterPatch, slavePatch;
        label dirIndex1, dirIndex2, dirIndex3, joinNumber;

        for (label glueI = 0; glueI < nglue; glueI++)
        {
            cfxFile >> masterPatch >> slavePatch;
            cfxFile >> dirIndex1 >> dirIndex2 >> dirIndex3 >> joinNumber;

            glueMasterPatches[glueI] = masterPatch - 1;
            glueSlavePatches[glueI] = slavePatch - 1;
        }
    }

    Info<< "Reading block points" << endl;

    forAll(blocks, blockI)
    {
        Info<< "block " << blockI << " is a ";
        blocks[blockI].readPoints(cfxFile);
    }

    Info<< "Calculating block offsets" << endl;

    labelList blockOffsets(nblock, -1);

    blockOffsets[0] = 0;

    label nMeshPoints = blocks[0].nBlockPoints();
    label nMeshCells = blocks[0].nBlockCells();

    for (label blockI = 1; blockI < nblock; blockI++)
    {
        nMeshPoints += blocks[blockI].nBlockPoints();
        nMeshCells +=  blocks[blockI].nBlockCells();

        blockOffsets[blockI] =
            blockOffsets[blockI - 1]
          + blocks[blockI - 1].nBlockPoints();
    }

    Info<< "Assembling patches" << endl;

    faceListList rawPatches(npatch);

    forAll(rawPatches, patchI)
    {
        const word& patchType = cfxPatchTypes[patchI];

        // reject volume patches
        if
        (
            patchType == "POROUS" || patchType == "SOLID"
         || patchType == "SOLCON" || patchType == "USER3D"
        )
        {
            patchMasterBlocks[patchI] = -1;
            rawPatches[patchI].setSize(0);
        }
        else
        {
            // read and create a 2-D patch
            rawPatches[patchI] =
                blocks[patchMasterBlocks[patchI]].patchFaces
                (
                    patchDirections[patchI],
                    patchRanges[patchI]
                );

        }
    }

    Info<< "Merging points ";

    labelList pointMergeList(nMeshPoints, -1);

    // In order to ensure robust merging, it is necessary to traverse
    // the patch glueing list until the pointMergeList stops changing.
    //

    // For efficiency, create merge pairs in the first pass
    labelListListList glueMergePairs(glueMasterPatches.size());

    forAll(glueMasterPatches, glueI)
    {
        const label masterPatch = glueMasterPatches[glueI];
        const label slavePatch = glueSlavePatches[glueI];

        const label blockPlabel = patchMasterBlocks[masterPatch];
        const label blockNlabel = patchMasterBlocks[slavePatch];

        const pointField& blockPpoints = blocks[blockPlabel].points();
        const pointField& blockNpoints = blocks[blockNlabel].points();

        const faceList& blockPFaces = rawPatches[masterPatch];
        const faceList& blockNFaces = rawPatches[slavePatch];

        labelListList& curPairs = glueMergePairs[glueI];
        curPairs.setSize(blockPFaces.size());

        if (blockPFaces.size() != blockNFaces.size())
        {
            FatalErrorIn(args.executable())
                << "Inconsistent number of faces for glue pair "
                << glueI << " between blocks " << blockPlabel + 1
                << " and " << blockNlabel + 1
                << abort(FatalError);
        }

        // Calculate sqr of the merge tolerance as 1/10th of the min
        // sqr point to point distance on the block face.  This is an
        // N^2 algorithm, sorry but I cannot quickly come up with
        // something better.

        scalar sqrMergeTol = GREAT;

        forAll(blockPFaces, blockPFaceLabel)
        {
            const labelList& blockPFacePoints =
                blockPFaces[blockPFaceLabel];

            forAll(blockPFacePoints, blockPFacePointI)
            {
                forAll(blockPFacePoints, blockPFacePointI2)
                {
                    if (blockPFacePointI != blockPFacePointI2)
                    {
                        sqrMergeTol =
                            min
                            (
                                sqrMergeTol,
                                magSqr
                                (
                                    blockPpoints
                                        [blockPFacePoints[blockPFacePointI]]
                                  - blockPpoints
                                        [blockPFacePoints[blockPFacePointI2]]
                                )
                            );
                    }
                }
            }
        }

        sqrMergeTol /= 10.0;

        register bool found = false;

        // N-squared point search over all points of all faces of
        // master block over all point of all faces of slave block
        forAll(blockPFaces, blockPFaceLabel)
        {
            const labelList& blockPFacePoints =
                blockPFaces[blockPFaceLabel];

            labelList& cp = curPairs[blockPFaceLabel];
            cp.setSize(blockPFacePoints.size());

        forAll(blockPFacePoints, blockPFacePointI)
        {
            found = false;

            forAll(blockNFaces, blockNFaceLabel)
            {
                const labelList& blockNFacePoints =
                    blockNFaces[blockNFaceLabel];

            forAll(blockNFacePoints, blockNFacePointI)
            {
                if
                (
                    magSqr
                    (
                        blockPpoints
                            [blockPFacePoints[blockPFacePointI]]
                      - blockNpoints
                            [blockNFacePoints[blockNFacePointI]]
                    )
                  < sqrMergeTol
                )
                {
                    // Found a new pair
                    found = true;

                    cp[blockPFacePointI] =
                        blockNFacePoints[blockNFacePointI];

                    label PpointLabel =
                        blockPFacePoints[blockPFacePointI]
                      + blockOffsets[blockPlabel];

                    label NpointLabel =
                        blockNFacePoints[blockNFacePointI]
                      + blockOffsets[blockNlabel];

                    label minPN = min(PpointLabel, NpointLabel);

                    if (pointMergeList[PpointLabel] != -1)
                    {
                        minPN = min(minPN, pointMergeList[PpointLabel]);
                    }

                    if (pointMergeList[NpointLabel] != -1)
                    {
                        minPN = min(minPN, pointMergeList[NpointLabel]);
                    }

                    pointMergeList[PpointLabel]
                  = pointMergeList[NpointLabel]
                  = minPN;

                    break;
                }
            }
            if (found) break;
            }
        }
        }
    }


    register bool changedPointMerge = false;
    label nPasses = 0;

    do
    {
        changedPointMerge = false;
        nPasses++;

        forAll(glueMasterPatches, glueI)
        {
            const label masterPatch = glueMasterPatches[glueI];
            const label slavePatch = glueSlavePatches[glueI];

            const label blockPlabel = patchMasterBlocks[masterPatch];
            const label blockNlabel = patchMasterBlocks[slavePatch];

            const faceList& blockPFaces = rawPatches[masterPatch];

            const labelListList& curPairs = glueMergePairs[glueI];

            forAll(blockPFaces, blockPFaceLabel)
            {
                const labelList& blockPFacePoints =
                    blockPFaces[blockPFaceLabel];

                const labelList& cp = curPairs[blockPFaceLabel];

                forAll(cp, blockPFacePointI)
                {
                    label PpointLabel =
                        blockPFacePoints[blockPFacePointI]
                      + blockOffsets[blockPlabel];

                    label NpointLabel =
                        cp[blockPFacePointI]
                      + blockOffsets[blockNlabel];

                    if
                    (
                        pointMergeList[PpointLabel]
                     != pointMergeList[NpointLabel]
                    )
                    {
                        changedPointMerge = true;

                        pointMergeList[PpointLabel]
                      = pointMergeList[NpointLabel]
                      = min
                        (
                            pointMergeList[PpointLabel],
                            pointMergeList[NpointLabel]
                        );
                    }
                }
            }
        }
        Info<< "." << flush;
    }
    while (changedPointMerge && nPasses < 8);
    Info<< endl;

    if (changedPointMerge == true)
    {
        FatalErrorIn(args.executable())
            << "Point merging failed after max number of passes."
            << abort(FatalError);
    }


    forAll(glueMasterPatches, glueI)
    {
        const label masterPatch = glueMasterPatches[glueI];
        const label slavePatch = glueSlavePatches[glueI];

        const label blockPlabel = patchMasterBlocks[masterPatch];
        const label blockNlabel = patchMasterBlocks[slavePatch];

        const faceList& blockPFaces = rawPatches[masterPatch];
        const faceList& blockNFaces = rawPatches[slavePatch];


        forAll(blockPFaces, blockPFaceLabel)
        {
            const labelList& blockPFacePoints
                = blockPFaces[blockPFaceLabel];

            forAll(blockPFacePoints, blockPFacePointI)
            {
                label PpointLabel =
                    blockPFacePoints[blockPFacePointI]
                  + blockOffsets[blockPlabel];

                if (pointMergeList[PpointLabel] == -1)
                {
                    FatalErrorIn(args.executable())
                        << "Unable to merge point " << blockPFacePointI
                        << " of face " << blockPFaceLabel
                        << " of block " << blockPlabel
                        << abort(FatalError);
                }
            }
        }

        forAll(blockNFaces, blockNFaceLabel)
        {
            const labelList& blockNFacePoints
                = blockNFaces[blockNFaceLabel];

            forAll(blockNFacePoints, blockNFacePointI)
            {
                label NpointLabel =
                    blockNFacePoints[blockNFacePointI]
                  + blockOffsets[blockNlabel];

                if (pointMergeList[NpointLabel] == -1)
                {
                    FatalErrorIn(args.executable())
                        << "Unable to merge point " << blockNFacePointI
                        << " of face " << blockNFaceLabel
                        << " of block " << blockNlabel
                        << abort(FatalError);
                }
            }
        }
    }


    // sort merge list to return new point label (in new shorter list)
    // given old point label
    label nNewPoints = 0;

    forAll(pointMergeList, pointLabel)
    {
        if (pointMergeList[pointLabel] > pointLabel)
        {
            FatalErrorIn(args.executable())
                << "ouch" << abort(FatalError);
        }

        if
        (
            (pointMergeList[pointLabel] == -1)
         || pointMergeList[pointLabel] == pointLabel
        )
        {
            pointMergeList[pointLabel] = nNewPoints;
            nNewPoints++;
        }
        else
        {
            pointMergeList[pointLabel] =
                pointMergeList[pointMergeList[pointLabel]];
        }
    }

    nMeshPoints = nNewPoints;

    Info<< "Creating points" << endl;

    pointField points(nMeshPoints);

    forAll(blocks, blockI)
    {
        const pointField& blockPoints = blocks[blockI].points();

        forAll(blockPoints, blockPointLabel)
        {
            points
            [
                pointMergeList
                [
                    blockPointLabel
                  + blockOffsets[blockI]
                ]
            ] = blockPoints[blockPointLabel];
        }
    }

    // Scale the points
    if (scaleFactor > 1.0 + SMALL || scaleFactor < 1.0 - SMALL)
    {
        points *= scaleFactor;
    }

    Info<< "Creating cells" << endl;

    cellShapeList cellShapes(nMeshCells);

    const cellModel& hex = *(cellModeller::lookup("hex"));

    label nCreatedCells = 0;

    forAll(blocks, blockI)
    {
        labelListList curBlockCells = blocks[blockI].blockCells();

        forAll(curBlockCells, blockCellI)
        {
            labelList cellPoints(curBlockCells[blockCellI].size());

            forAll(cellPoints, pointI)
            {
                cellPoints[pointI] =
                    pointMergeList
                    [
                        curBlockCells[blockCellI][pointI]
                      + blockOffsets[blockI]
                    ];
            }

            cellShapes[nCreatedCells] = cellShape(hex, cellPoints);

            nCreatedCells++;
        }
    }

    Info<< "Creating boundary patches" << endl;

    faceListList boundary(npatch);
    wordList patchNames(npatch);
    wordList patchTypes(npatch);
    word defaultFacesName = "defaultFaces";
    word defaultFacesType = wallPolyPatch::typeName;

    label nCreatedPatches = 0;

    forAll(rawPatches, patchI)
    {
        if (rawPatches[patchI].size() && cfxPatchTypes[patchI] != "BLKBDY")
        {
            // Check if this name has been already created
            label existingPatch = -1;

            for (label oldPatchI = 0; oldPatchI < nCreatedPatches; oldPatchI++)
            {
                if (patchNames[oldPatchI] == cfxPatchNames[patchI])
                {
                    existingPatch = oldPatchI;
                    break;
                }
            }

            const faceList& curRawPatch = rawPatches[patchI];
            label curBlock = patchMasterBlocks[patchI];

            if (existingPatch >= 0)
            {
                Info<< "CFX patch " << patchI
                    << ", of type " << cfxPatchTypes[patchI]
                    << ", name " << cfxPatchNames[patchI]
                    << " already exists as OpenFOAM patch " << existingPatch
                    << ".  Adding faces." << endl;

                faceList& renumberedPatch = boundary[existingPatch];
                label oldSize = renumberedPatch.size();
                renumberedPatch.setSize(oldSize + curRawPatch.size());

                forAll(curRawPatch, faceI)
                {
                    const face& oldFace = curRawPatch[faceI];

                    face& newFace = renumberedPatch[oldSize + faceI];
                    newFace.setSize(oldFace.size());

                    forAll(oldFace, pointI)
                    {
                        newFace[pointI] =
                            pointMergeList
                            [
                                oldFace[pointI]
                              + blockOffsets[curBlock]
                            ];
                    }
                }
            }
            else
            {
                // Real patch to be created
                faceList& renumberedPatch = boundary[nCreatedPatches];
                renumberedPatch.setSize(curRawPatch.size());

                forAll(curRawPatch, faceI)
                {
                    const face& oldFace = curRawPatch[faceI];

                    face& newFace = renumberedPatch[faceI];
                    newFace.setSize(oldFace.size());

                    forAll(oldFace, pointI)
                    {
                        newFace[pointI] =
                            pointMergeList
                            [
                                oldFace[pointI]
                              + blockOffsets[curBlock]
                            ];
                    }
                }

                Info<< "CFX patch " << patchI
                    << ", of type " << cfxPatchTypes[patchI]
                    << ", name " << cfxPatchNames[patchI]
                    << " converted into OpenFOAM patch " << nCreatedPatches
                    << " type ";

                if (cfxPatchTypes[patchI] == "WALL")
                {
                    Info<< "wall." << endl;

                    patchTypes[nCreatedPatches] = wallPolyPatch::typeName;
                    patchNames[nCreatedPatches] = cfxPatchNames[patchI];
                    nCreatedPatches++;
                }
                else if (cfxPatchTypes[patchI] == "SYMMET")
                {
                    Info<< "symmetryPlane." << endl;

                    patchTypes[nCreatedPatches] = symmetryPolyPatch::typeName;
                    patchNames[nCreatedPatches] = cfxPatchNames[patchI];
                    nCreatedPatches++;
                }
                else if
                (
                    cfxPatchTypes[patchI] == "INLET"
                 || cfxPatchTypes[patchI] == "OUTLET"
                 || cfxPatchTypes[patchI] == "PRESS"
                 || cfxPatchTypes[patchI] == "CNDBDY"
                 || cfxPatchTypes[patchI] == "USER2D"
                )
                {
                    Info<< "generic." << endl;

                    patchTypes[nCreatedPatches] = polyPatch::typeName;
                    patchNames[nCreatedPatches] = cfxPatchNames[patchI];
                    nCreatedPatches++;
                }
                else
                {
                    FatalErrorIn(args.executable())
                        << "Unrecognised CFX patch type "
                        << cfxPatchTypes[patchI]
                        << abort(FatalError);
                }
            }
        }
    }

    boundary.setSize(nCreatedPatches);
    patchTypes.setSize(nCreatedPatches);
    patchNames.setSize(nCreatedPatches);

    PtrList<dictionary> patchDicts;

    preservePatchTypes
    (
        runTime,
        runTime.constant(),
        polyMesh::meshSubDir,
        patchNames,
        patchDicts,
        defaultFacesName,
        defaultFacesType
    );

    // Add information to dictionary
    forAll(patchNames, patchI)
    {
        if (!patchDicts.set(patchI))
        {
            patchDicts.set(patchI, new dictionary());
        }
        // Add but not overwrite
        patchDicts[patchI].add("type", patchTypes[patchI], false);
    }

    polyMesh pShapeMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        xferMove(points),
        cellShapes,
        boundary,
        patchNames,
        patchDicts,
        defaultFacesName,
        defaultFacesType
    );

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    Info<< "Writing polyMesh" << endl;
    pShapeMesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
