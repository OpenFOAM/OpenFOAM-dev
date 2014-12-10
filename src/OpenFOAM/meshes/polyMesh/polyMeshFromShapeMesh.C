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

Description
    Create polyMesh from cell and patch shapes

\*---------------------------------------------------------------------------*/

#include "polyMesh.H"
#include "Time.H"
#include "primitiveMesh.H"
#include "DynamicList.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelListList Foam::polyMesh::cellShapePointCells
(
    const cellShapeList& c
) const
{
    List<DynamicList<label, primitiveMesh::cellsPerPoint_> >
        pc(points().size());

    // For each cell
    forAll(c, i)
    {
        // For each vertex
        const labelList& labels = c[i];

        forAll(labels, j)
        {
            // Set working point label
            label curPoint = labels[j];
            DynamicList<label, primitiveMesh::cellsPerPoint_>& curPointCells =
                pc[curPoint];

            // Enter the cell label in the point's cell list
            curPointCells.append(i);
        }
    }

    labelListList pointCellAddr(pc.size());

    forAll(pc, pointI)
    {
        pointCellAddr[pointI].transfer(pc[pointI]);
    }

    return pointCellAddr;
}


Foam::labelList Foam::polyMesh::facePatchFaceCells
(
    const faceList& patchFaces,
    const labelListList& pointCells,
    const faceListList& cellsFaceShapes,
    const label patchID
) const
{
    register bool found;

    labelList FaceCells(patchFaces.size());

    forAll(patchFaces, fI)
    {
        found = false;

        const face& curFace = patchFaces[fI];
        const labelList& facePoints = patchFaces[fI];

        forAll(facePoints, pointI)
        {
            const labelList& facePointCells = pointCells[facePoints[pointI]];

            forAll(facePointCells, cellI)
            {
                faceList cellFaces = cellsFaceShapes[facePointCells[cellI]];

                forAll(cellFaces, cellFace)
                {
                    if (cellFaces[cellFace] == curFace)
                    {
                        // Found the cell corresponding to this face
                        FaceCells[fI] = facePointCells[cellI];

                        found = true;
                    }
                    if (found) break;
                }
                if (found) break;
            }
            if (found) break;
        }

        if (!found)
        {
            FatalErrorIn
            (
                "polyMesh::facePatchFaceCells(const faceList& patchFaces,"
                "const labelListList& pointCells,"
                "const faceListList& cellsFaceShapes,"
                "const label patchID)"
            )   << "face " << fI << " in patch " << patchID
                << " does not have neighbour cell"
                << " face: " << patchFaces[fI]
                << abort(FatalError);
        }
    }

    return FaceCells;
}


//- Set faces_, calculate cells and patchStarts.
void Foam::polyMesh::setTopology
(
    const cellShapeList& cellsAsShapes,
    const faceListList& boundaryFaces,
    const wordList& boundaryPatchNames,
    labelList& patchSizes,
    labelList& patchStarts,
    label& defaultPatchStart,
    label& nFaces,
    cellList& cells
)
{
    // Calculate the faces of all cells
    // Initialise maximum possible numer of mesh faces to 0
    label maxFaces = 0;

    // Set up a list of face shapes for each cell
    faceListList cellsFaceShapes(cellsAsShapes.size());
    cells.setSize(cellsAsShapes.size());

    forAll(cellsFaceShapes, cellI)
    {
        cellsFaceShapes[cellI] = cellsAsShapes[cellI].faces();

        cells[cellI].setSize(cellsFaceShapes[cellI].size());

        // Initialise cells to -1 to flag undefined faces
        static_cast<labelList&>(cells[cellI]) = -1;

        // Count maximum possible numer of mesh faces
        maxFaces += cellsFaceShapes[cellI].size();
    }

    // Set size of faces array to maximum possible number of mesh faces
    faces_.setSize(maxFaces);

    // Initialise number of faces to 0
    nFaces = 0;

    // set reference to point-cell addressing
    labelListList PointCells = cellShapePointCells(cellsAsShapes);

    bool found = false;

    forAll(cells, cellI)
    {
        // Note:
        // Insertion cannot be done in one go as the faces need to be
        // added into the list in the increasing order of neighbour
        // cells.  Therefore, all neighbours will be detected first
        // and then added in the correct order.

        const faceList& curFaces = cellsFaceShapes[cellI];

        // Record the neighbour cell
        labelList neiCells(curFaces.size(), -1);

        // Record the face of neighbour cell
        labelList faceOfNeiCell(curFaces.size(), -1);

        label nNeighbours = 0;

        // For all faces ...
        forAll(curFaces, faceI)
        {
            // Skip faces that have already been matched
            if (cells[cellI][faceI] >= 0) continue;

            found = false;

            const face& curFace = curFaces[faceI];

            // Get the list of labels
            const labelList& curPoints = curFace;

            // For all points
            forAll(curPoints, pointI)
            {
                // dGget the list of cells sharing this point
                const labelList& curNeighbours =
                    PointCells[curPoints[pointI]];

                // For all neighbours
                forAll(curNeighbours, neiI)
                {
                    label curNei = curNeighbours[neiI];

                    // Reject neighbours with the lower label
                    if (curNei > cellI)
                    {
                        // Get the list of search faces
                        const faceList& searchFaces = cellsFaceShapes[curNei];

                        forAll(searchFaces, neiFaceI)
                        {
                            if (searchFaces[neiFaceI] == curFace)
                            {
                                // Match!!
                                found = true;

                                // Record the neighbour cell and face
                                neiCells[faceI] = curNei;
                                faceOfNeiCell[faceI] = neiFaceI;
                                nNeighbours++;

                                break;
                            }
                        }
                        if (found) break;
                    }
                    if (found) break;
                }
                if (found) break;
            } // End of current points
        }  // End of current faces

        // Add the faces in the increasing order of neighbours
        for (label neiSearch = 0; neiSearch < nNeighbours; neiSearch++)
        {
            // Find the lowest neighbour which is still valid
            label nextNei = -1;
            label minNei = cells.size();

            forAll(neiCells, ncI)
            {
                if (neiCells[ncI] > -1 && neiCells[ncI] < minNei)
                {
                    nextNei = ncI;
                    minNei = neiCells[ncI];
                }
            }

            if (nextNei > -1)
            {
                // Add the face to the list of faces
                faces_[nFaces] = curFaces[nextNei];

                // Set cell-face and cell-neighbour-face to current face label
                cells[cellI][nextNei] = nFaces;
                cells[neiCells[nextNei]][faceOfNeiCell[nextNei]] = nFaces;

                // Stop the neighbour from being used again
                neiCells[nextNei] = -1;

                // Increment number of faces counter
                nFaces++;
            }
            else
            {
                FatalErrorIn
                (
                    "polyMesh::setTopology\n"
                    "(\n"
                    "    const cellShapeList& cellsAsShapes,\n"
                    "    const faceListList& boundaryFaces,\n"
                    "    const wordList& boundaryPatchNames,\n"
                    "    labelList& patchSizes,\n"
                    "    labelList& patchStarts,\n"
                    "    label& defaultPatchStart,\n"
                    "    label& nFaces,\n"
                    "    cellList& cells\n"
                    ")"
                )   << "Error in internal face insertion"
                    << abort(FatalError);
            }
        }
    }

    // Do boundary faces

    patchSizes.setSize(boundaryFaces.size(), -1);
    patchStarts.setSize(boundaryFaces.size(), -1);

    forAll(boundaryFaces, patchI)
    {
        const faceList& patchFaces = boundaryFaces[patchI];

        labelList curPatchFaceCells =
            facePatchFaceCells
            (
                patchFaces,
                PointCells,
                cellsFaceShapes,
                patchI
            );

        // Grab the start label
        label curPatchStart = nFaces;

        forAll(patchFaces, faceI)
        {
            const face& curFace = patchFaces[faceI];

            const label cellInside = curPatchFaceCells[faceI];

            faces_[nFaces] = curFace;

            // get faces of the cell inside
            const faceList& facesOfCellInside = cellsFaceShapes[cellInside];

            bool found = false;

            forAll(facesOfCellInside, cellFaceI)
            {
                if (facesOfCellInside[cellFaceI] == curFace)
                {
                    if (cells[cellInside][cellFaceI] >= 0)
                    {
                        FatalErrorIn
                        (
                            "polyMesh::setTopology\n"
                            "(\n"
                            "    const cellShapeList& cellsAsShapes,\n"
                            "    const faceListList& boundaryFaces,\n"
                            "    const wordList& boundaryPatchNames,\n"
                            "    labelList& patchSizes,\n"
                            "    labelList& patchStarts,\n"
                            "    label& defaultPatchStart,\n"
                            "    label& nFaces,\n"
                            "    cellList& cells\n"
                            ")"
                        )   << "Trying to specify a boundary face " << curFace
                            << " on the face on cell " << cellInside
                            << " which is either an internal face or already "
                            << "belongs to some other patch.  This is face "
                            << faceI << " of patch "
                            << patchI << " named "
                            << boundaryPatchNames[patchI] << "."
                            << abort(FatalError);
                    }

                    found = true;

                    cells[cellInside][cellFaceI] = nFaces;

                    break;
                }
            }

            if (!found)
            {
                FatalErrorIn("polyMesh::polyMesh(... construct from shapes...)")
                    << "face " << faceI << " of patch " << patchI
                    << " does not seem to belong to cell " << cellInside
                    << " which, according to the addressing, "
                    << "should be next to it."
                    << abort(FatalError);
            }

            // increment the counter of faces
            nFaces++;
        }

        patchSizes[patchI] = nFaces - curPatchStart;
        patchStarts[patchI] = curPatchStart;
    }

    // Grab "non-existing" faces and put them into a default patch

    defaultPatchStart = nFaces;

    forAll(cells, cellI)
    {
        labelList& curCellFaces = cells[cellI];

        forAll(curCellFaces, faceI)
        {
            if (curCellFaces[faceI] == -1) // "non-existent" face
            {
                curCellFaces[faceI] = nFaces;
                faces_[nFaces] = cellsFaceShapes[cellI][faceI];

                nFaces++;
            }
        }
    }

    // Reset the size of the face list
    faces_.setSize(nFaces);

    return ;
}


Foam::polyMesh::polyMesh
(
    const IOobject& io,
    const Xfer<pointField>& points,
    const cellShapeList& cellsAsShapes,
    const faceListList& boundaryFaces,
    const wordList& boundaryPatchNames,
    const wordList& boundaryPatchTypes,
    const word& defaultBoundaryPatchName,
    const word& defaultBoundaryPatchType,
    const wordList& boundaryPatchPhysicalTypes,
    const bool syncPar
)
:
    objectRegistry(io),
    primitiveMesh(),
    points_
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        points
    ),
    faces_
    (
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    owner_
    (
        IOobject
        (
            "owner",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    neighbour_
    (
        IOobject
        (
            "neighbour",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    clearedPrimitives_(false),
    boundary_
    (
        IOobject
        (
            "boundary",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        boundaryFaces.size() + 1    // add room for a default patch
    ),
    bounds_(points_, syncPar),
    comm_(UPstream::worldComm),
    geometricD_(Vector<label>::zero),
    solutionD_(Vector<label>::zero),
    tetBasePtIsPtr_(NULL),
    cellTreePtr_(NULL),
    pointZones_
    (
        IOobject
        (
            "pointZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    faceZones_
    (
        IOobject
        (
            "faceZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    cellZones_
    (
        IOobject
        (
            "cellZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    globalMeshDataPtr_(NULL),
    moving_(false),
    topoChanging_(false),
    curMotionTimeIndex_(time().timeIndex()),
    oldPointsPtr_(NULL)
{
    if (debug)
    {
        Info<<"Constructing polyMesh from cell and boundary shapes." << endl;
    }

    // Remove all of the old mesh files if they exist
    removeFiles(instance());

    // Calculate faces and cells
    labelList patchSizes;
    labelList patchStarts;
    label defaultPatchStart;
    label nFaces;
    cellList cells;
    setTopology
    (
        cellsAsShapes,
        boundaryFaces,
        boundaryPatchNames,
        patchSizes,
        patchStarts,
        defaultPatchStart,
        nFaces,
        cells
    );

    // Warning: Patches can only be added once the face list is
    // completed, as they hold a subList of the face list
    forAll(boundaryFaces, patchI)
    {
        // add the patch to the list
        boundary_.set
        (
            patchI,
            polyPatch::New
            (
                boundaryPatchTypes[patchI],
                boundaryPatchNames[patchI],
                patchSizes[patchI],
                patchStarts[patchI],
                patchI,
                boundary_
            )
        );

        if
        (
            boundaryPatchPhysicalTypes.size()
         && boundaryPatchPhysicalTypes[patchI].size()
        )
        {
            boundary_[patchI].physicalType() =
                boundaryPatchPhysicalTypes[patchI];
        }
    }

    label nAllPatches = boundaryFaces.size();


    label nDefaultFaces = nFaces - defaultPatchStart;
    if (syncPar)
    {
        reduce(nDefaultFaces, sumOp<label>());
    }

    if (nDefaultFaces > 0)
    {
        WarningIn("polyMesh::polyMesh(... construct from shapes...)")
            << "Found " << nDefaultFaces
            << " undefined faces in mesh; adding to default patch." << endl;

        // Check if there already exists a defaultFaces patch as last patch
        // and reuse it.
        label patchI = findIndex(boundaryPatchNames, defaultBoundaryPatchName);

        if (patchI != -1)
        {
            if (patchI != boundaryFaces.size()-1 || boundary_[patchI].size())
            {
                FatalErrorIn("polyMesh::polyMesh(... construct from shapes...)")
                    << "Default patch " << boundary_[patchI].name()
                    << " already has faces in it or is not"
                    << " last in list of patches." << exit(FatalError);
            }

            WarningIn("polyMesh::polyMesh(... construct from shapes...)")
                << "Reusing existing patch " << patchI
                << " for undefined faces." << endl;

            boundary_.set
            (
                patchI,
                polyPatch::New
                (
                    boundaryPatchTypes[patchI],
                    boundaryPatchNames[patchI],
                    nFaces - defaultPatchStart,
                    defaultPatchStart,
                    patchI,
                    boundary_
                )
            );
        }
        else
        {
            boundary_.set
            (
                nAllPatches,
                polyPatch::New
                (
                    defaultBoundaryPatchType,
                    defaultBoundaryPatchName,
                    nFaces - defaultPatchStart,
                    defaultPatchStart,
                    boundary_.size() - 1,
                    boundary_
                )
            );

            nAllPatches++;
        }
    }

    // Reset the size of the boundary
    boundary_.setSize(nAllPatches);

    // Set the primitive mesh
    initMesh(cells);

    if (syncPar)
    {
        // Calculate topology for the patches (processor-processor comms etc.)
        boundary_.updateMesh();

        // Calculate the geometry for the patches (transformation tensors etc.)
        boundary_.calcGeometry();
    }

    if (debug)
    {
        if (checkMesh())
        {
            Info<< "Mesh OK" << endl;
        }
    }
}


Foam::polyMesh::polyMesh
(
    const IOobject& io,
    const Xfer<pointField>& points,
    const cellShapeList& cellsAsShapes,
    const faceListList& boundaryFaces,
    const wordList& boundaryPatchNames,
    const PtrList<dictionary>& boundaryDicts,
    const word& defaultBoundaryPatchName,
    const word& defaultBoundaryPatchType,
    const bool syncPar
)
:
    objectRegistry(io),
    primitiveMesh(),
    points_
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        points
    ),
    faces_
    (
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    owner_
    (
        IOobject
        (
            "owner",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    neighbour_
    (
        IOobject
        (
            "neighbour",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    clearedPrimitives_(false),
    boundary_
    (
        IOobject
        (
            "boundary",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        boundaryFaces.size() + 1    // add room for a default patch
    ),
    bounds_(points_, syncPar),
    comm_(UPstream::worldComm),
    geometricD_(Vector<label>::zero),
    solutionD_(Vector<label>::zero),
    tetBasePtIsPtr_(NULL),
    cellTreePtr_(NULL),
    pointZones_
    (
        IOobject
        (
            "pointZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    faceZones_
    (
        IOobject
        (
            "faceZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    cellZones_
    (
        IOobject
        (
            "cellZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    globalMeshDataPtr_(NULL),
    moving_(false),
    topoChanging_(false),
    curMotionTimeIndex_(time().timeIndex()),
    oldPointsPtr_(NULL)
{
    if (debug)
    {
        Info<<"Constructing polyMesh from cell and boundary shapes." << endl;
    }

    // Remove all of the old mesh files if they exist
    removeFiles(instance());

    // Calculate faces and cells
    labelList patchSizes;
    labelList patchStarts;
    label defaultPatchStart;
    label nFaces;
    cellList cells;
    setTopology
    (
        cellsAsShapes,
        boundaryFaces,
        boundaryPatchNames,
        patchSizes,
        patchStarts,
        defaultPatchStart,
        nFaces,
        cells
    );

    // Warning: Patches can only be added once the face list is
    // completed, as they hold a subList of the face list
    forAll(boundaryDicts, patchI)
    {
        dictionary patchDict(boundaryDicts[patchI]);

        patchDict.set("nFaces", patchSizes[patchI]);
        patchDict.set("startFace", patchStarts[patchI]);

        // add the patch to the list
        boundary_.set
        (
            patchI,
            polyPatch::New
            (
                boundaryPatchNames[patchI],
                patchDict,
                patchI,
                boundary_
            )
        );
    }

    label nAllPatches = boundaryFaces.size();

    label nDefaultFaces = nFaces - defaultPatchStart;
    if (syncPar)
    {
        reduce(nDefaultFaces, sumOp<label>());
    }

    if (nDefaultFaces > 0)
    {
        WarningIn("polyMesh::polyMesh(... construct from shapes...)")
            << "Found " << nDefaultFaces
            << " undefined faces in mesh; adding to default patch." << endl;

        // Check if there already exists a defaultFaces patch as last patch
        // and reuse it.
        label patchI = findIndex(boundaryPatchNames, defaultBoundaryPatchName);

        if (patchI != -1)
        {
            if (patchI != boundaryFaces.size()-1 || boundary_[patchI].size())
            {
                FatalErrorIn("polyMesh::polyMesh(... construct from shapes...)")
                    << "Default patch " << boundary_[patchI].name()
                    << " already has faces in it or is not"
                    << " last in list of patches." << exit(FatalError);
            }

            WarningIn("polyMesh::polyMesh(... construct from shapes...)")
                << "Reusing existing patch " << patchI
                << " for undefined faces." << endl;

            boundary_.set
            (
                patchI,
                polyPatch::New
                (
                    boundary_[patchI].type(),
                    boundary_[patchI].name(),
                    nFaces - defaultPatchStart,
                    defaultPatchStart,
                    patchI,
                    boundary_
                )
            );
        }
        else
        {
            boundary_.set
            (
                nAllPatches,
                polyPatch::New
                (
                    defaultBoundaryPatchType,
                    defaultBoundaryPatchName,
                    nFaces - defaultPatchStart,
                    defaultPatchStart,
                    boundary_.size() - 1,
                    boundary_
                )
            );

            nAllPatches++;
        }
    }

    // Reset the size of the boundary
    boundary_.setSize(nAllPatches);

    // Set the primitive mesh
    initMesh(cells);

    if (syncPar)
    {
        // Calculate topology for the patches (processor-processor comms etc.)
        boundary_.updateMesh();

        // Calculate the geometry for the patches (transformation tensors etc.)
        boundary_.calcGeometry();
    }

    if (debug)
    {
        if (checkMesh())
        {
            Info << "Mesh OK" << endl;
        }
    }
}


// ************************************************************************* //
