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
    find the tetrahedron in which the position is.
    while searching for the tet, store the tet
    closest to the position.
    This way, if position is not inside any tet,
    choose the closest one.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool interpolationCellPointFace<Type>::findTet
(
    const vector& position,
    const label nFace,
    vector tetPoints[],
    label tetLabelCandidate[],
    label tetPointLabels[],
    scalar phi[],
    scalar phiCandidate[],
    label& closestFace,
    scalar& minDistance
) const
{
    bool foundTet = false;

    const labelList& thisFacePoints = this->pMeshFaces_[nFace];
    tetPoints[2] = this->pMeshFaceCentres_[nFace];

    label pointi = 0;

    while (pointi < thisFacePoints.size() && !foundTet)
    {
        label nextPointLabel = (pointi + 1) % thisFacePoints.size();

        tetPointLabels[0] = thisFacePoints[pointi];
        tetPointLabels[1] = thisFacePoints[nextPointLabel];

        tetPoints[0] = this->pMeshPoints_[tetPointLabels[0]];
        tetPoints[1] = this->pMeshPoints_[tetPointLabels[1]];

        bool inside = true;
        scalar dist = 0.0;

        for (label n=0; n<4; n++)
        {
            label p1 = (n + 1) % 4;
            label p2 = (n + 2) % 4;
            label p3 = (n + 3) % 4;

            vector referencePoint, faceNormal;
            referencePoint = tetPoints[p1];

            faceNormal =
                (tetPoints[p3] - tetPoints[p1])
              ^ (tetPoints[p2] - tetPoints[p1]);

            faceNormal /= mag(faceNormal);

            // correct normal to point into the tet
            vector v0 = tetPoints[n] - referencePoint;
            scalar correct = v0 & faceNormal;
            if (correct < 0)
            {
                faceNormal = -faceNormal;
            }

            vector v1 = position - referencePoint + SMALL*faceNormal;
            scalar rightSide = v1 & faceNormal;

            // since normal is inwards, inside the tet
            // is defined as all dot-products being positive
            inside = inside && (rightSide >= 0);

            scalar phiLength = (position - referencePoint) & faceNormal;

            scalar maxLength =
                max(VSMALL, (tetPoints[n] - referencePoint) & faceNormal);

            phi[n] = phiLength/maxLength;

            dist += phi[n];
        }

        if (!inside)
        {
            if (mag(dist - 1.0) < minDistance)
            {
                minDistance = mag(dist - 1.0);
                closestFace = nFace;

                for (label i=0; i<4; i++)
                {
                    phiCandidate[i] = phi[i];
                }

                tetLabelCandidate[0] = tetPointLabels[0];
                tetLabelCandidate[1] = tetPointLabels[1];
            }
        }

        foundTet = inside;

        pointi++;
    }

    if (foundTet)
    {
        closestFace = nFace;
    }

    return foundTet;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
