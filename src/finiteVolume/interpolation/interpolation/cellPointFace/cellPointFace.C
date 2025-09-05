/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "cellPointFace.H"
#include "linear.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::interpolations::cellPointFace<Type>::findTet
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

    const labelList& thisFacePoints = this->mesh_.faces()[nFace];
    tetPoints[2] = this->mesh_.faceCentres()[nFace];

    label pointi = 0;

    while (pointi < thisFacePoints.size() && !foundTet)
    {
        label nextPointLabel = (pointi + 1) % thisFacePoints.size();

        tetPointLabels[0] = thisFacePoints[pointi];
        tetPointLabels[1] = thisFacePoints[nextPointLabel];

        tetPoints[0] = this->mesh_.points()[tetPointLabels[0]];
        tetPoints[1] = this->mesh_.points()[tetPointLabels[1]];

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

            vector v1 = position - referencePoint + small*faceNormal;
            scalar rightSide = v1 & faceNormal;

            // since normal is inwards, inside the tet
            // is defined as all dot-products being positive
            inside = inside && (rightSide >= 0);

            scalar phiLength = (position - referencePoint) & faceNormal;

            scalar maxLength =
                max(vSmall, (tetPoints[n] - referencePoint) & faceNormal);

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


template<class Type>
bool Foam::interpolations::cellPointFace<Type>::findTriangle
(
    const vector& position,
    const label nFace,
    label tetPointLabels[],
    scalar phi[]
) const
{
    bool foundTriangle = false;
    vector tetPoints[3];
    const labelList& facePoints = this->mesh_.faces()[nFace];
    tetPoints[2] = this->mesh_.faceCentres()[nFace];

    label pointi = 0;

    while (pointi < facePoints.size() && !foundTriangle)
    {
        // The triangle is constructed from face center and one
        // face edge
        label nextPointLabel = (pointi + 1) % facePoints.size();

        tetPointLabels[0] = facePoints[pointi];
        tetPointLabels[1] = facePoints[nextPointLabel];

        tetPoints[0] = this->mesh_.points()[tetPointLabels[0]];
        tetPoints[1] = this->mesh_.points()[tetPointLabels[1]];

        vector fc = (tetPoints[0] + tetPoints[1] + tetPoints[2])/3.0;

        vector newPos = position + small*(fc-position);

        // calculate triangle edge vectors and triangle face normal
        // the 'i':th edge is opposite node i
        vector edge[3], normal[3];
        for (label i=0; i<3; i++)
        {
            label ip0 = (i+1) % 3;
            label ipp = (i+2) % 3;
            edge[i] = tetPoints[ipp]-tetPoints[ip0];
        }

        vector triangleFaceNormal = edge[1] ^ edge[2];

        // calculate edge normal (pointing inwards)
        for (label i=0; i<3; i++)
        {
            normal[i] = triangleFaceNormal ^ edge[i];
            normal[i] /= mag(normal[i]) + vSmall;
        }

        // check if position is inside triangle
        bool inside = true;
        for (label i=0; i<3; i++)
        {
            label ip = (i+1) % 3;
            inside = inside && (((newPos - tetPoints[ip]) & edge[i]) >= 0);
        }

        if (inside)
        {
            foundTriangle = true;

            // calculate phi-values
            for (label i=0; i<3; i++)
            {
                label ip = (i+1) % 3;
                scalar phiMax = max(vSmall, normal[i] & edge[ip]);
                scalar phiLength = (position-tetPoints[ip]) & normal[i];
                phi[i] = phiLength/phiMax;
            }
        }

        pointi++;
    }

    return foundTriangle;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolations::cellPointFace<Type>::cellPointFace
(
    const VolField<Type>& psi
)
:
    fieldInterpolation<Type, cellPointFace<Type>>(psi),
    volPointInterpolation<Type>(psi),
    psis_(linearInterpolate(psi))
{}


template<class Type>
Foam::interpolations::cellPointFace<Type>::cellPointFace
(
    const cellPointFace<Type>& i
)
:
    fieldInterpolation<Type, cellPointFace<Type>>(i),
    volPointInterpolation<Type>(i),
    psis_(i.psis_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::interpolations::cellPointFace<Type>::interpolate
(
    const vector& position,
    const label celli,
    const label facei
) const
{
    Type ts[4];
    vector tetPoints[4];
    scalar phi[4], phiCandidate[4];
    label tetLabelCandidate[2], tetPointLabels[2];

    Type t = Zero;

    // only use face information when the position is on a face
    if (facei < 0)
    {
        const vector& cellCentre = this->mesh_.cellCentres()[celli];
        const labelList& cellFaces = this->mesh_.cells()[celli];

        vector projection = position - cellCentre;
        tetPoints[3] = cellCentre;

        // *********************************************************************
        // project the cell-center through the point onto the face
        // and get the closest face, ...
        // *********************************************************************

        bool foundTet = false;
        label closestFace = -1;
        scalar minDistance = great;

        forAll(cellFaces, facei)
        {
            label nFace = cellFaces[facei];

            vector normal = this->mesh_.faceAreas()[nFace];
            normal /= mag(normal);

            const vector& faceCentreTmp = this->mesh_.faceCentres()[nFace];

            scalar multiplierNumerator = (faceCentreTmp - cellCentre) & normal;
            scalar multiplierDenominator = projection & normal;

            // if normal and projection are not orthogonal this could
            // be the one...
            if (mag(multiplierDenominator) > small)
            {
                scalar multiplier = multiplierNumerator/multiplierDenominator;
                vector iPoint = cellCentre + multiplier*projection;
                scalar dist = mag(position - iPoint);

                if (dist < minDistance)
                {
                    closestFace = nFace;
                    minDistance = dist;
                }
            }
        }

        // *********************************************************************
        // find the tetrahedron containing 'position'
        // from the cell center, face center and
        // two other points on the face
        // *********************************************************************

        minDistance = great;
        if (closestFace != -1)
        {
            label nFace = closestFace;
            foundTet = findTet
            (
                position,
                nFace,
                tetPoints,
                tetLabelCandidate,
                tetPointLabels,
                phi,
                phiCandidate,
                closestFace,
                minDistance
            );
        }

        if (!foundTet)
        {
            // check if the position is 'just' outside a tet
            if (minDistance < 1.0e-5)
            {
                foundTet = true;
                for (label i=0; i<4; i++)
                {
                    phi[i] = phiCandidate[i];
                }
                tetPointLabels[0] = tetLabelCandidate[0];
                tetPointLabels[1] = tetLabelCandidate[1];
            }
        }

        // *********************************************************************
        // if the search failed check all the cell-faces
        // *********************************************************************

        if (!foundTet)
        {
            minDistance = great;

            label facei = 0;
            while (facei < cellFaces.size() && !foundTet)
            {
                label nFace = cellFaces[facei];
                if (nFace < this->mesh_.faceAreas().size())
                {
                    foundTet = findTet
                    (
                        position,
                        nFace,
                        tetPoints,
                        tetLabelCandidate,
                        tetPointLabels,
                        phi,
                        phiCandidate,
                        closestFace,
                        minDistance
                    );
                }
                facei++;
            }
        }

        if (!foundTet)
        {
            // check if the position is 'just' outside a tet
            // this time with a more tolerant limit
            if (minDistance < 1.0e-3)
            {
                foundTet = true;
                for (label i=0; i<4; i++)
                {
                    phi[i] = phiCandidate[i];
                }
                tetPointLabels[0] = tetLabelCandidate[0];
                tetPointLabels[1] = tetLabelCandidate[1];
            }
        }

        // *********************************************************************
        // if the tet was found do the interpolation,
        // otherwise use the closest face value
        // *********************************************************************

        if (foundTet)
        {
            for (label i=0; i<2; i++)
            {
                ts[i] = this->psip_[tetPointLabels[i]];
            }

            if (closestFace < psis_.size())
            {
                ts[2] = psis_[closestFace];
            }
            else
            {
                label patchi =
                    this->mesh_.boundaryMesh().whichPatch(closestFace);

                // If the boundary patch is not empty use the face value
                // else use the cell value
                if (this->psi_.boundaryField()[patchi].size())
                {
                    ts[2] = this->psi_.boundaryField()[patchi]
                    [
                        this->mesh_.boundaryMesh()[patchi].whichFace
                        (
                            closestFace
                        )
                    ];
                }
                else
                {
                    ts[2] = this->psi_[celli];
                }
            }

            ts[3] = this->psi_[celli];

            for (label n=0; n<4; n++)
            {
                phi[n] = min(1.0, phi[n]);
                phi[n] = max(0.0, phi[n]);

                t += phi[n]*ts[n];
            }
        }
        else
        {
            InfoInFunction
                << "search failed; using closest cellFace value" << endl
                << "cell number " << celli << tab
                << "position " << position << endl;

            if (closestFace < psis_.size())
            {
                t = psis_[closestFace];
            }
            else
            {
                label patchi =
                    this->mesh_.boundaryMesh().whichPatch(closestFace);

                // If the boundary patch is not empty use the face value
                // else use the cell value
                if (this->psi_.boundaryField()[patchi].size())
                {
                    t = this->psi_.boundaryField()[patchi]
                    [
                        this->mesh_.boundaryMesh()[patchi].whichFace
                        (
                            closestFace
                        )
                    ];
                }
                else
                {
                    t = this->psi_[celli];
                }
            }
        }
    }
    else
    {
        bool foundTriangle = findTriangle
        (
            position,
            facei,
            tetPointLabels,
            phi
        );

        if (foundTriangle)
        {
            // add up the point values ...
            for (label i=0; i<2; i++)
            {
                Type vel = this->psip_[tetPointLabels[i]];
                t += phi[i]*vel;
            }

            // ... and the face value
            if (facei < psis_.size())
            {
                t += phi[2]*psis_[facei];
            }
            else
            {
                label patchi = this->mesh_.boundaryMesh().whichPatch(facei);

                // If the boundary patch is not empty use the face value
                // else use the cell value
                if (this->psi_.boundaryField()[patchi].size())
                {
                    t += phi[2]*this->psi_.boundaryField()[patchi]
                        [this->mesh_.boundaryMesh()[patchi].whichFace(facei)];
                }
                else
                {
                    t += phi[2]*this->psi_[celli];
                }
            }
        }
        else
        {
            // use face value only
            if (facei < psis_.size())
            {
                t = psis_[facei];
            }
            else
            {
                label patchi = this->mesh_.boundaryMesh().whichPatch(facei);

                // If the boundary patch is not empty use the face value
                // else use the cell value
                if (this->psi_.boundaryField()[patchi].size())
                {
                    t = this->psi_.boundaryField()[patchi]
                        [this->mesh_.boundaryMesh()[patchi].whichFace(facei)];
                }
                else
                {
                    t = this->psi_[celli];
                }
            }
        }
    }

    return t;
}


// ************************************************************************* //
