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
    Construct an extruded triangular prism cell shape from three straight edges

\*---------------------------------------------------------------------------*/

#include "cellShapeRecognition.H"
#include "labelList.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

cellShape extrudedTriangleCellShape
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
    static const cellModel* prismModelPtr_ = NULL;

    if (!prismModelPtr_)
    {
        prismModelPtr_ = cellModeller::lookup("prism");
    }

    const cellModel& prism = *prismModelPtr_;

    // Checking
    if (faceLabels.size() != 3)
    {
        FatalErrorIn
        (
            "extrudedTriangleCellShape(const label cellIndex, "
            "const labelList& faceLabels, const faceList& faces, "
            "const labelList& owner, const labelList& neighbour, "
            "const label pointOffset, faceList& frontAndBackFaces)"
        )   << "Trying to create a triangle with " << faceLabels.size()
            << " faces"
            << abort(FatalError);
    }

    // make a list of outward-pointing faces
    labelListList localFaces(3);

    forAll(faceLabels, faceI)
    {
        const label curFaceLabel = faceLabels[faceI];

        const face& curFace = faces[curFaceLabel];

        if (curFace.size() != 2)
        {
            FatalErrorIn
            (
                "extrudedTriangleCellShape(const label cellIndex, "
                "const labelList& faceLabels, const faceList& faces, "
                "const labelList& owner, const labelList& neighbour, "
                "const label pointOffset, faceList& frontAndBackFaces)"
            )   << "face " << curFaceLabel
                << "does not have 2 vertices. Number of vertices: " << curFace
                << abort(FatalError);
        }

        if (owner[curFaceLabel] == cellIndex)
        {
            localFaces[faceI] = curFace;
        }
        else if (neighbour[curFaceLabel] == cellIndex)
        {
            // Reverse the face.  Note: it is necessary to reverse by
            // hand to preserve connectivity of a 2-D mesh.
            //
            localFaces[faceI].setSize(curFace.size());

            forAllReverse(curFace, i)
            {
                localFaces[faceI][curFace.size() - i - 1] =
                    curFace[i];
            }
        }
        else
        {
            FatalErrorIn
            (
                "extrudedTriangleCellShape(const label cellIndex, "
                "const labelList& faceLabels, const faceList& faces, "
                "const labelList& owner, const labelList& neighbour, "
                "const label pointOffset, faceList& frontAndBackFaces)"
            )   << "face " << curFaceLabel
                << " does not belong to cell " << cellIndex
                << ". Face owner: " << owner[curFaceLabel] << " neighbour: "
                << neighbour[curFaceLabel]
                << abort(FatalError);
        }
    }

    // Create a label list for the model
    if (localFaces[0][1] == localFaces[1][0])
    {
        // Set front and back plane faces
        labelList missingPlaneFace(3);

        // front plane
        missingPlaneFace[0] = localFaces[0][0];
        missingPlaneFace[1] = localFaces[1][1];
        missingPlaneFace[2] = localFaces[0][1];

        frontAndBackFaces[2*cellIndex] = face(missingPlaneFace);

        // back plane
        missingPlaneFace[0] = localFaces[0][0] + pointOffset;
        missingPlaneFace[1] = localFaces[0][1] + pointOffset;
        missingPlaneFace[2] = localFaces[1][1] + pointOffset;

        frontAndBackFaces[2*cellIndex + 1] = face(missingPlaneFace);

        // make a cell
        labelList cellShapeLabels(6);

        cellShapeLabels[0] = localFaces[0][0];
        cellShapeLabels[1] = localFaces[0][1];
        cellShapeLabels[2] = localFaces[1][1];

        cellShapeLabels[3] = localFaces[0][0] + pointOffset;
        cellShapeLabels[4] = localFaces[0][1] + pointOffset;
        cellShapeLabels[5] = localFaces[1][1] + pointOffset;

        return cellShape(prism, cellShapeLabels);
    }
    else if (localFaces[0][1] == localFaces[2][0])
    {
        // Set front and back plane faces
        labelList missingPlaneFace(3);

        // front plane
        missingPlaneFace[0] = localFaces[0][0];
        missingPlaneFace[1] = localFaces[2][1];
        missingPlaneFace[2] = localFaces[0][1];

        frontAndBackFaces[2*cellIndex] = face(missingPlaneFace);

        // back plane
        missingPlaneFace[0] = localFaces[0][0] + pointOffset;
        missingPlaneFace[1] = localFaces[0][1] + pointOffset;
        missingPlaneFace[2] = localFaces[2][1] + pointOffset;

        frontAndBackFaces[2*cellIndex + 1] = face(missingPlaneFace);

        // make a cell
        labelList cellShapeLabels(6);

        cellShapeLabels[0] = localFaces[0][0];
        cellShapeLabels[1] = localFaces[0][1];
        cellShapeLabels[2] = localFaces[2][1];

        cellShapeLabels[3] = localFaces[0][0] + pointOffset;
        cellShapeLabels[4] = localFaces[0][1] + pointOffset;
        cellShapeLabels[5] = localFaces[2][1] + pointOffset;

        return cellShape(prism, cellShapeLabels);
    }
    else
    {
        FatalErrorIn
        (
            "extrudedTriangleCellShape(const label cellIndex, "
            "const labelList& faceLabels, const faceList& faces, "
            "const labelList& owner, const labelList& neighbour, "
            "const label pointOffset, faceList& frontAndBackFaces)"
        )   << "Problem with edge matching. Edges: " << localFaces
            << abort(FatalError);
    }

    // Return added to keep compiler happy
    return cellShape(prism, labelList(0));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
