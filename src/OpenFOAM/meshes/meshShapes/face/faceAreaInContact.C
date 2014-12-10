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

\*---------------------------------------------------------------------------*/

#include "face.H"
#include "scalarField.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Calculate area in contact given displacement of vertices relative to
// the face plane. Positive displacement is above the face (no contact);
// negative is in contact
Foam::scalar Foam::face::areaInContact
(
    const pointField& meshPoints,
    const scalarField& v
) const
{
    // Assemble the vertex values
    const labelList& labels = *this;

    scalarField vertexValue(labels.size());

    forAll(labels, i)
    {
        vertexValue[i] = v[labels[i]];
    }


    // Loop through vertexValue. If all greater that 0 return 0 (no contact);
    // if all less than zero return 1
    // all zeros is assumed to be in contact.

    bool allPositive = true;
    bool allNegative = true;

    forAll(vertexValue, vI)
    {
        if (vertexValue[vI] > 0)
        {
            allNegative = false;
        }
        else
        {
            allPositive = false;
        }
    }

    if (allPositive)
    {
        return 0.0;
    }

    if (allNegative)
    {
        return 1.0;
    }

    // There is a partial contact.
    // Algorithm:
    // Go through all edges. if both vertex values for the edge are
    // positive, discard. If one is positive and one is negative,
    // create a point and start the edge with it. If both are
    // negative, add the edge into the new face.  When finished,
    // calculate area of new face and return relative area (0<x<1)

    // Dimension new point list to max possible size
    const labelList& faceLabels = *this;

    pointField newFacePoints(2*size());
    label nNewFacePoints = 0;

    for (label vI = 0; vI < size() - 1; vI++)
    {
        if (vertexValue[vI] <= 0)
        {
            // This is a point in contact
            newFacePoints[nNewFacePoints] = meshPoints[faceLabels[vI]];
            nNewFacePoints++;
        }

        if
        (
            (vertexValue[vI] > 0 && vertexValue[vI + 1] < 0)
         || (vertexValue[vI] < 0 && vertexValue[vI + 1] > 0)
        )
        {
            // Edge intersection. Calculate intersection point and add to list
            point intersection =
                meshPoints[faceLabels[vI]]
              + vertexValue[vI]/(vertexValue[vI + 1] - vertexValue[vI])
                *(meshPoints[faceLabels[vI]] - meshPoints[faceLabels[vI + 1]]);

            newFacePoints[nNewFacePoints] = intersection;
            nNewFacePoints++;
        }
    }

    // Do last point by hand
    if (vertexValue[size() - 1] <= 0)
    {
        // This is a point in contact
        newFacePoints[nNewFacePoints] = meshPoints[faceLabels[size() - 1]];
        nNewFacePoints++;
    }

    if
    (
        (vertexValue[size() - 1] > 0 && vertexValue[0] < 0)
     || (vertexValue[size() - 1] < 0 && vertexValue[0] > 0)
    )
    {
        // Edge intersection. Calculate intersection point and add to list
        point intersection =
            meshPoints[faceLabels[size() - 1]]
          + vertexValue[size() - 1]/(vertexValue[0] - vertexValue[size() - 1])
            *(meshPoints[faceLabels[size() - 1]] - meshPoints[faceLabels[0]]);

        newFacePoints[nNewFacePoints] = intersection;
        nNewFacePoints++;
    }

    newFacePoints.setSize(nNewFacePoints);

    // Make a labelList for the sub-face (points are ordered!)
    labelList sfl(newFacePoints.size());

    forAll(sfl, sflI)
    {
        sfl[sflI] = sflI;
    }

    // Calculate relative area
    return face(sfl).mag(newFacePoints)/(mag(meshPoints) + VSMALL);
}


// ************************************************************************* //
