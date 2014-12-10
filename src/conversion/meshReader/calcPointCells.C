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
    calculate point cells - ie, the cells attached to each point

    - remove unused points, adjust pointCells and cellFaces accordingly
\*---------------------------------------------------------------------------*/

#include "meshReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::meshReader::calcPointCells() const
{
    static const label UNIT_POINT_CELLS = 12;

    if (pointCellsPtr_)
    {
        FatalErrorIn("meshReader::calcPointCells() const")
            << "pointCells already calculated"
            << abort(FatalError);
    }

    label nPoints = points_.size();

    pointCellsPtr_ = new labelListList(nPoints);
    labelListList& ptCells = *pointCellsPtr_;

    forAll(ptCells, i)
    {
        ptCells[i].setSize(UNIT_POINT_CELLS);
    }

    // Initialize the list of labels which will hold the count of the
    // actual number of cells per point during the analysis
    labelList cellCount(nPoints, 0);

    // Note. Unlike the standard point-cell algorithm, which asks the cell for
    // the supporting point labels, we need to work based on the cell faces.
    // This is because some of the faces do not come from the cell shape.
    // It is also advantageous to remove duplicates from the point-cell
    // addressing, because this removes a lot of waste later.

    faceListList& cFaces = cellFaces();

    // For each cell
    forAll(cFaces, cellI)
    {
        const faceList& faces = cFaces[cellI];

        forAll(faces, i)
        {
            // For each vertex
            const labelList& labels = faces[i];

            forAll(labels, j)
            {
                // Set working point label
                label curPoint = labels[j];
                labelList& curPointCells = ptCells[curPoint];
                label curCount = cellCount[curPoint];

                // check if the cell has been added before
                bool found = false;

                for (label f = 0; f < curCount; f++)
                {
                    if (curPointCells[f] == cellI)
                    {
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    // If the list of pointCells is not big enough, double it
                    if (curPointCells.size() <= curCount)
                    {
                        curPointCells.setSize(curPointCells.size()*2);
                    }

                    // Enter the cell label in the point's cell list
                    curPointCells[curCount] = cellI;

                    // Increment the cell count for the point addressed
                    cellCount[curPoint]++;
                }
            }
        }
    }

    // report and remove unused points
    // - adjust points, pointCells, and cellFaces accordingly
    label pointI = 0;
    labelList oldToNew(nPoints, -1);

    forAll(ptCells, i)
    {
        ptCells[i].setSize(cellCount[i]);
        if (cellCount[i] > 0)
        {
            oldToNew[i] = pointI++;
        }
    }

    // report unused points
    if (nPoints > pointI)
    {
        Info<< "removing " << (nPoints - pointI) << " unused points" << endl;

        nPoints = pointI;

        // adjust points and truncate - bend const-ness
        pointField& adjustedPoints = const_cast<pointField&>(points_);

        inplaceReorder(oldToNew, adjustedPoints);
        adjustedPoints.setSize(nPoints);

        // adjust pointCells and truncate
        inplaceReorder(oldToNew, ptCells);
        ptCells.setSize(nPoints);

        // adjust cellFaces - this could be faster
        // For each cell
        forAll(cFaces, cellI)
        {
            faceList& faces = cFaces[cellI];

            // For each face
            forAll(faces, i)
            {
                inplaceRenumber(oldToNew, faces[i]);
            }
        }
    }
}


const Foam::labelListList& Foam::meshReader::pointCells() const
{
    if (!pointCellsPtr_)
    {
        calcPointCells();
    }

    return *pointCellsPtr_;
}


// ************************************************************************* //
