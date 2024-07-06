/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Construct an extruded hex cell shape from four straight edges

\*---------------------------------------------------------------------------*/

#include "cellShapeRecognition.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

cellShape extrudedQuadCellShape
(
    const label cellIndex,
    const labelList& faceLabels,
    const faceList& faces,
    const labelList& owner,
    const labelList& neighbour,
    const label pointOffset,
    faceList& frontAndBackFaces
)
{
    static const cellModel* hexModelPtr_ = nullptr;

    if (!hexModelPtr_)
    {
        hexModelPtr_ = cellModeller::lookup("hex");
    }

    const cellModel& hex = *hexModelPtr_;

    // Checking
    if (faceLabels.size() != 4)
    {
        FatalErrorInFunction
            << "Trying to create a quad with " << faceLabels.size() << " faces"
            << abort(FatalError);
    }

    // make a list of outward-pointing faces
    labelListList localFaces(4);

    forAll(faceLabels, facei)
    {
        const label curFaceLabel = faceLabels[facei];

        const face& curFace = faces[curFaceLabel];

        if (curFace.size() != 2)
        {
            FatalErrorInFunction
                << "face " << curFaceLabel
                << "does not have 2 vertices. Number of vertices: " << curFace
                << abort(FatalError);
        }

        if (owner[curFaceLabel] == cellIndex)
        {
            localFaces[facei] = curFace;
        }
        else if (neighbour[curFaceLabel] == cellIndex)
        {
            // Reverse the face.  Note: it is necessary to reverse by
            // hand to preserve connectivity of a 2-D mesh.
            //
            localFaces[facei].setSize(curFace.size());

            forAllReverse(curFace, i)
            {
                localFaces[facei][curFace.size() - i - 1] =
                    curFace[i];
            }
        }
        else
        {
            FatalErrorInFunction
                << "face " << curFaceLabel
                << " does not belong to cell " << cellIndex
                << ". Face owner: " << owner[curFaceLabel] << " neighbour: "
                << neighbour[curFaceLabel]
                << abort(FatalError);
        }
    }

    // Create a label list for the model
    // This is done by finding two edges that do not share any vertices.
    // Knowing the opposite pair of edges (with normals pointing outward
    // is enough to make a cell
    if
    (
        localFaces[0][0] != localFaces[1][1]
     && localFaces[0][1] != localFaces[1][0]
    )
    {
        // Set front and back plane faces
        labelList missingPlaneFace(4);

        // front plane
        missingPlaneFace[0] = localFaces[0][0];
        missingPlaneFace[1] = localFaces[1][1];
        missingPlaneFace[2] = localFaces[1][0];
        missingPlaneFace[3] = localFaces[0][1];

        frontAndBackFaces[2*cellIndex] = face(missingPlaneFace);

        // back plane
        missingPlaneFace[0] = localFaces[0][0] + pointOffset;
        missingPlaneFace[1] = localFaces[0][1] + pointOffset;
        missingPlaneFace[2] = localFaces[1][0] + pointOffset;
        missingPlaneFace[3] = localFaces[1][1] + pointOffset;

        frontAndBackFaces[2*cellIndex + 1] = face(missingPlaneFace);

        // make a cell
        labelList cellShapeLabels(8);

        cellShapeLabels[0] = localFaces[0][0];
        cellShapeLabels[1] = localFaces[0][1];
        cellShapeLabels[2] = localFaces[1][0];
        cellShapeLabels[3] = localFaces[1][1];

        cellShapeLabels[4] = localFaces[0][0] + pointOffset;
        cellShapeLabels[5] = localFaces[0][1] + pointOffset;
        cellShapeLabels[6] = localFaces[1][0] + pointOffset;
        cellShapeLabels[7] = localFaces[1][1] + pointOffset;


        return cellShape(hex, cellShapeLabels);
    }
    else if
    (
        localFaces[0][0] != localFaces[2][1]
     && localFaces[0][1] != localFaces[2][0]
    )
    {
        // Set front and back plane faces
        labelList missingPlaneFace(4);

        // front plane
        missingPlaneFace[0] = localFaces[0][0];
        missingPlaneFace[1] = localFaces[2][1];
        missingPlaneFace[2] = localFaces[2][0];
        missingPlaneFace[3] = localFaces[0][1];

        frontAndBackFaces[2*cellIndex] = face(missingPlaneFace);

        // back plane
        missingPlaneFace[0] = localFaces[0][0] + pointOffset;
        missingPlaneFace[1] = localFaces[0][1] + pointOffset;
        missingPlaneFace[2] = localFaces[2][0] + pointOffset;
        missingPlaneFace[3] = localFaces[2][1] + pointOffset;

        frontAndBackFaces[2*cellIndex + 1] = face(missingPlaneFace);

        // make a cell
        labelList cellShapeLabels(8);

        cellShapeLabels[0] = localFaces[0][0];
        cellShapeLabels[1] = localFaces[0][1];
        cellShapeLabels[2] = localFaces[2][0];
        cellShapeLabels[3] = localFaces[2][1];

        cellShapeLabels[4] = localFaces[0][0] + pointOffset;
        cellShapeLabels[5] = localFaces[0][1] + pointOffset;
        cellShapeLabels[6] = localFaces[2][0] + pointOffset;
        cellShapeLabels[7] = localFaces[2][1] + pointOffset;


        return cellShape(hex, cellShapeLabels);
    }
    else if
    (
        localFaces[0][0] != localFaces[3][1]
     && localFaces[0][1] != localFaces[3][0]
    )
    {
        // Set front and back plane faces
        labelList missingPlaneFace(4);

        // front plane
        missingPlaneFace[0] = localFaces[0][0];
        missingPlaneFace[1] = localFaces[3][1];
        missingPlaneFace[2] = localFaces[3][0];
        missingPlaneFace[3] = localFaces[0][1];

        frontAndBackFaces[2*cellIndex] = face(missingPlaneFace);

        // back plane
        missingPlaneFace[0] = localFaces[0][0] + pointOffset;
        missingPlaneFace[1] = localFaces[0][1] + pointOffset;
        missingPlaneFace[2] = localFaces[3][0] + pointOffset;
        missingPlaneFace[3] = localFaces[3][1] + pointOffset;

        frontAndBackFaces[2*cellIndex + 1] = face(missingPlaneFace);

        // make a cell
        labelList cellShapeLabels(8);

        cellShapeLabels[0] = localFaces[0][0];
        cellShapeLabels[1] = localFaces[0][1];
        cellShapeLabels[2] = localFaces[3][0];
        cellShapeLabels[3] = localFaces[3][1];

        cellShapeLabels[4] = localFaces[0][0] + pointOffset;
        cellShapeLabels[5] = localFaces[0][1] + pointOffset;
        cellShapeLabels[6] = localFaces[3][0] + pointOffset;
        cellShapeLabels[7] = localFaces[3][1] + pointOffset;


        return cellShape(hex, cellShapeLabels);
    }
    else
    {
        FatalErrorInFunction
            << "Problem with edge matching. Edges: " << localFaces
            << abort(FatalError);
    }

    // Return added to keep compiler happy
    return cellShape(hex, labelList(0));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
