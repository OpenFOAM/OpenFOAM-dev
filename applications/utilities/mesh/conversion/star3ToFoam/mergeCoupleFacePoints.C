/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Merge duplicate points created by arbitrary face coupling and remove unused
    points.

\*---------------------------------------------------------------------------*/

#include "starMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void starMesh::mergeCoupleFacePoints()
{
    // mark all used points by looping through all faces in two goes.
    // First, go into every cell and find min edge length.  Use a
    // fraction of that as a merge tolerance.  Loop through all the
    // points of the cell and query a merge against every other point
    // of the same cell based on the current tolerance.  If two points
    // merge, find out if any of them already belongs to a "merge set"
    // (group of points that merge together.  If so, add the current
    // points into the merge set and make sure that the merge label
    // for the whole set is the lowest of all encountered vertices.
    // This is stored in a mergeSetID list, under the index of the
    // merge set.  Once all cells (and thus points) are visited, go
    // through the renumbering list and for each merging point use the
    // label of the merge set as the new point label.
    // This is VERY fancy. Use care if/when changing.

    Info<< endl << "Creating merge sets" << endl;

    // Create a renumbering list for points
    // In the first instance the renumbering list is used as a
    // mergeSetID storage
    labelList renumberPoints(points_.size(), -1);

    // create storage of mergeSetIDs
    label mergeIncrement = points_.size()/10;
    labelList mergeSetID(mergeIncrement, -1);

    label nMergeSets = 0;

    forAll(cellFaces_, cellI)
    {
        const faceList& curFaces = cellFaces_[cellI];

        // create a list of all points for the cell with duplicates
        // and find the shortest edge length

        label nPointsInCell = 0;

        scalar pointMergeTol = GREAT;

        forAll(curFaces, faceI)
        {
            nPointsInCell += curFaces[faceI].size();

            edgeList curEdges = curFaces[faceI].edges();

            forAll(curEdges, edgeI)
            {
                scalar length = curEdges[edgeI].mag(points_);

                if (length < pointMergeTol)
                {
                    pointMergeTol = length;
                }
            }
        }

        // set merge tolerance as a fraction of the shortest edge
        pointMergeTol /= 100.0;

        // create a list of points
        labelList cellPoints(nPointsInCell);
        label nAddedPoints = 0;

        forAll(curFaces, faceI)
        {
            const face& f = curFaces[faceI];

            forAll(f, fI)
            {
                cellPoints[nAddedPoints] = f[fI];
                nAddedPoints++;
            }
        }

        // loop n-squared through the local points and merge them up
        for (label firstPI = 0; firstPI < cellPoints.size() - 1; firstPI++)
        {
            for
            (
                label otherPI = firstPI;
                otherPI < cellPoints.size();
                otherPI++
            )
            {
                if (cellPoints[firstPI] != cellPoints[otherPI])
                {
                    label a = cellPoints[firstPI];
                    label b = cellPoints[otherPI];

                    if (edge (a, b).mag(points_) < pointMergeTol)
                    {
                        // found a pair of points to merge
#                       ifdef DEBUG_MERGE
                        Info<< "Merging points " << a << " and " << b << endl;
#                       endif

                        // are the two points in a merge group?
                        label mergeSetA = -1;
                        label mergeSetB = -1;

                        if (renumberPoints[a] > -1)
                        {
                            mergeSetA = renumberPoints[a];
                        }

                        if (renumberPoints[b] > -1)
                        {
                            mergeSetB = renumberPoints[b];
                        }

                        if (mergeSetA == -1 && mergeSetB == -1)
                        {
                            // add new merge group
#                           ifdef DEBUG_MERGE
                            Info<< "adding new merge group " << nMergeSets
                                << endl;
#                           endif

                            // mark points as belonging to a new merge set
                            renumberPoints[a] = nMergeSets;
                            renumberPoints[b] = nMergeSets;

                            mergeSetID[nMergeSets] = min(a, b);
                            nMergeSets++;

                            if (nMergeSets >= mergeSetID.size())
                            {
                                Info<< "Resizing mergeSetID" << endl;

                                mergeSetID.setSize
                                    (mergeSetID.size() + mergeIncrement);
                            }
                        }
                        else if (mergeSetA == -1 && mergeSetB != -1)
                        {
#                           ifdef DEBUG_MERGE
                            Info<< "adding point a into the merge set of b. "
                                << "a: " << a << endl;
#                           endif

                            // add point a into the merge set of b
                            renumberPoints[a] = mergeSetB;

                            // reset the min label of the mergeSet
                            mergeSetID[mergeSetB] =
                                min(a, mergeSetID[mergeSetB]);
                        }
                        else if  (mergeSetA != -1 && mergeSetB == -1)
                        {
#                           ifdef DEBUG_MERGE
                            Info<< "adding point b into the merge set of a. "
                                << "b: " << b << endl;
#                           endif

                            // add point b into the merge set of a
                            renumberPoints[b] = mergeSetA;

                            // reset the min label of the mergeSet
                            mergeSetID[mergeSetA] =
                                min(b, mergeSetID[mergeSetA]);
                        }
                        else if (mergeSetA != mergeSetB)
                        {
                            // Points already belong to two different merge
                            // sets. Eliminate the higher merge set
                            label minMerge = min(mergeSetA, mergeSetB);
                            label maxMerge = max(mergeSetA, mergeSetB);

#                           ifdef DEBUG_MERGE
                            Info<< "Points already belong to two "
                                << "different merge sets. "
                                << "Eliminate the higher merge set. Sets: "
                                << minMerge << " and " << maxMerge << endl;
#                           endif

                            forAll(renumberPoints, elimI)
                            {
                                if (renumberPoints[elimI] == maxMerge)
                                {
                                    renumberPoints[elimI] = minMerge;
                                }
                            }

                            // set the lower label
                            mergeSetID[minMerge] =
                                min(mergeSetID[minMerge], mergeSetID[maxMerge]);
                        }
                    }
                }
            }
        }
    }

    mergeSetID.setSize(nMergeSets);

    Info<< "Finished creating merge sets.  Number of merge sets: "
        << nMergeSets << "." << endl;

    // Insert the primary point renumbering into the list
    // Take care of possibly unused points in the list
    forAll(renumberPoints, pointI)
    {
        if (renumberPoints[pointI] < 0)
        {
            // point not merged
            renumberPoints[pointI] = pointI;
        }
        else
        {
            renumberPoints[pointI] = mergeSetID[renumberPoints[pointI]];
        }
    }

    // Now every point carries either its own label or the lowest of
    // the labels of the points it is merged with. Now do PRELIMINARY
    // renumbering of all faces.  This will only be used to see which
    // points are still used!

    forAll(cellFaces_, cellI)
    {
        faceList& prelimFaces = cellFaces_[cellI];

        forAll(prelimFaces, faceI)
        {
            face oldFacePoints = prelimFaces[faceI];

            face& prelimFacePoints = prelimFaces[faceI];

            forAll(prelimFacePoints, pointI)
            {
                if (renumberPoints[oldFacePoints[pointI]] < 0)
                {
                    FatalErrorIn("starMesh::mergeCoupleFacePoints()")
                        << "Error in point renumbering. Old face: "
                        << oldFacePoints << endl
                        << "prelim face: " << prelimFacePoints
                        << abort(FatalError);
                }

                prelimFacePoints[pointI] =
                    renumberPoints[oldFacePoints[pointI]];
            }
        }
    }

    // First step complete. Reset the renumbering list, mark all used points,
    // re-create the point list and renumber the whole lot
    renumberPoints = 0;

    forAll(cellFaces_, cellI)
    {
        const faceList& curFaces = cellFaces_[cellI];

        forAll(curFaces, faceI)
        {
            const face& curFacePoints = curFaces[faceI];

            forAll(curFacePoints, pointI)
            {
                renumberPoints[curFacePoints[pointI]]++;
            }
        }
    }

    forAll(cellShapes_, cellI)
    {
        const labelList& curLabels = cellShapes_[cellI];

        forAll(curLabels, pointI)
        {
            if (renumberPoints[curLabels[pointI]] == 0)
            {
                FatalErrorIn("starMesh::mergeCoupleFacePoints()")
                    << "Error in point merging for cell "
                    << cellI << ". STAR index: " << starCellID_[cellI]
                    << ". " << endl
                    << "Point index: " << curLabels[pointI] << " STAR index "
                    << starPointID_[curLabels[pointI]] << endl
                    << "Please check the geometry for the cell." << endl;
            }
        }
    }

    label nUsedPoints = 0;

    forAll(renumberPoints, pointI)
    {
        if (renumberPoints[pointI] > 0)
        {
            // This should be OK as the compressed points list will always
            // have less points that the original lists.  Even if there is
            // no points removed, this will copy the list back onto itself
            //
            renumberPoints[pointI] = nUsedPoints;
            points_[nUsedPoints] = points_[pointI];

            nUsedPoints++;
        }
        else
        {
            renumberPoints[pointI] = -1;
        }
    }

    Info<< "Total number of points: " << points_.size() << endl
        << "Number of used points: " << nUsedPoints << endl;

    // reset number of points which need to be sorted
    points_.setSize(nUsedPoints);

    Info<< "Renumbering all faces" << endl;

    forAll(cellFaces_, cellI)
    {
        faceList& newFaces = cellFaces_[cellI];

        forAll(newFaces, faceI)
        {
            face oldFacePoints = newFaces[faceI];

            face& newFacePoints = newFaces[faceI];

            forAll(newFacePoints, pointI)
            {
                if (renumberPoints[oldFacePoints[pointI]] < 0)
                {
                    FatalErrorIn("starMesh::mergeCoupleFacePoints()")
                        << "Error in point renumbering for point "
                        << oldFacePoints[pointI]
                        << ". Renumbering index is -1." << endl
                        << "Old face: " << oldFacePoints << endl
                        << "New face: " << newFacePoints << abort(FatalError);
                }

                newFacePoints[pointI] = renumberPoints[oldFacePoints[pointI]];
            }
        }
    }

    Info<< "Renumbering all cell shapes" << endl;

    forAll(cellShapes_, cellI)
    {
        labelList oldLabels = cellShapes_[cellI];

        labelList& curLabels = cellShapes_[cellI];

        forAll(curLabels, pointI)
        {
            if (renumberPoints[curLabels[pointI]] < 0)
            {
                FatalErrorIn("starMesh::mergeCoupleFacePoints()")
                    << "Error in point renumbering for cell "
                    << cellI << ". STAR index: " << starCellID_[cellI]
                    << ". " << endl
                    << "Point index: " << curLabels[pointI] << " STAR index "
                    << starPointID_[curLabels[pointI]] << " returns invalid "
                    << "renumbering index: "
                    << renumberPoints[curLabels[pointI]] << "." << endl
                    << "Old cellShape:  " << oldLabels << endl
                    << "New cell shape: " << curLabels << abort(FatalError);
                }

            curLabels[pointI] = renumberPoints[oldLabels[pointI]];
        }
    }

    Info<< "Renumbering STAR point lookup" << endl;

    labelList oldStarPointID = starPointID_;

    starPointID_ = -1;

    forAll(starPointID_, pointI)
    {
        if (renumberPoints[pointI] > -1)
        {
            starPointID_[renumberPoints[pointI]] = oldStarPointID[pointI];
        }
    }

    // point-cell addressing has changed. Force it to be re-created
    deleteDemandDrivenData(pointCellsPtr_);
}


// ************************************************************************* //
