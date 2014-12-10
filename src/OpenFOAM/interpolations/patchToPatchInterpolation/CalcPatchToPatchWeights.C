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

#include "PatchToPatchInterpolation.H"
#include "objectHit.H"
#include "pointHit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class FromPatch, class ToPatch>
scalar PatchToPatchInterpolation<FromPatch, ToPatch>::projectionTol_ = 0.05;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FromPatch, class ToPatch>
void PatchToPatchInterpolation<FromPatch, ToPatch>::calcPointAddressing() const
{
    // Calculate pointWeights

    pointWeightsPtr_ = new FieldField<Field, scalar>(toPatch_.nPoints());
    FieldField<Field, scalar>& pointWeights = *pointWeightsPtr_;

    pointDistancePtr_ = new scalarField(toPatch_.nPoints(), GREAT);
    scalarField& pointDistance = *pointDistancePtr_;

    const pointField& fromPatchPoints = fromPatch_.localPoints();
    const List<typename FromPatch::FaceType>& fromPatchFaces =
        fromPatch_.localFaces();

    const pointField& toPatchPoints = toPatch_.localPoints();
    const vectorField& projectionDirection = toPatch_.pointNormals();
    const edgeList& toPatchEdges = toPatch_.edges();
    const labelListList& toPatchPointEdges = toPatch_.pointEdges();

    if (debug)
    {
        Info<< "projecting points" << endl;
    }

    List<objectHit> proj =
        toPatch_.projectPoints(fromPatch_, projectionDirection, alg_, dir_);

    pointAddressingPtr_ = new labelList(proj.size(), -1);
    labelList& pointAddressing = *pointAddressingPtr_;

    bool doWeights = false;

    forAll(pointAddressing, pointI)
    {
        doWeights = false;

        const typename FromPatch::FaceType& hitFace =
            fromPatchFaces[proj[pointI].hitObject()];

        point hitPoint = point::zero;

        if (proj[pointI].hit())
        {
            // A hit exists
            doWeights = true;

            pointAddressing[pointI] = proj[pointI].hitObject();

            pointHit curHit =
                hitFace.ray
                (
                    toPatchPoints[pointI],
                    projectionDirection[pointI],
                    fromPatchPoints,
                    alg_,
                    dir_
                );

            // Grab distance to target
            if (dir_ == intersection::CONTACT_SPHERE)
            {
                pointDistance[pointI] =
                    hitFace.contactSphereDiameter
                    (
                        toPatchPoints[pointI],
                        projectionDirection[pointI],
                        fromPatchPoints
                    );
            }
            else
            {
                pointDistance[pointI] = curHit.distance();
            }

            // Grab hit point
            hitPoint = curHit.hitPoint();
        }
        else if (projectionTol_ > SMALL)
        {
            // Check for a near miss
            pointHit ph =
                hitFace.ray
                (
                    toPatchPoints[pointI],
                    projectionDirection[pointI],
                    fromPatchPoints,
                    alg_,
                    dir_
                );

            scalar dist =
                Foam::mag
                (
                    toPatchPoints[pointI]
                  + projectionDirection[pointI]*ph.distance()
                  - ph.missPoint()
                );

            // Calculate the local tolerance
            scalar minEdgeLength = GREAT;

            // Do shortest edge of hit object
            edgeList hitFaceEdges =
                fromPatchFaces[proj[pointI].hitObject()].edges();

            forAll(hitFaceEdges, edgeI)
            {
                minEdgeLength =
                    min
                    (
                        minEdgeLength,
                        hitFaceEdges[edgeI].mag(fromPatchPoints)
                    );
            }

            const labelList& curEdges = toPatchPointEdges[pointI];

            forAll(curEdges, edgeI)
            {
                minEdgeLength =
                    min
                    (
                        minEdgeLength,
                        toPatchEdges[curEdges[edgeI]].mag(toPatchPoints)
                    );
            }

            if (dist < minEdgeLength*projectionTol_)
            {
                // This point is being corrected
                doWeights = true;

                pointAddressing[pointI] = proj[pointI].hitObject();

                // Grab nearest point on face as hit point
                hitPoint = ph.missPoint();

                // Grab distance to target
                if (dir_ == intersection::CONTACT_SPHERE)
                {
                    pointDistance[pointI] =
                        hitFace.contactSphereDiameter
                        (
                            toPatchPoints[pointI],
                            projectionDirection[pointI],
                            fromPatchPoints
                        );
                }
                else
                {
                    pointDistance[pointI] =
                        (
                            projectionDirection[pointI]
                            /mag(projectionDirection[pointI])
                        )
                      & (hitPoint - toPatchPoints[pointI]);
                }
            }
        }

        if (doWeights)
        {
            // Set interpolation pointWeights
            pointWeights.set(pointI, new scalarField(hitFace.size()));

            pointField hitFacePoints = hitFace.points(fromPatchPoints);

            forAll(hitFacePoints, masterPointI)
            {
                pointWeights[pointI][masterPointI] =
                    1.0/
                    (
                        mag
                        (
                            hitFacePoints[masterPointI]
                          - hitPoint
                        )
                      + VSMALL
                    );
            }

            pointWeights[pointI] /= sum(pointWeights[pointI]);
        }
        else
        {
            pointWeights.set(pointI, new scalarField(0));
        }
    }
}


template<class FromPatch, class ToPatch>
void PatchToPatchInterpolation<FromPatch, ToPatch>::calcFaceAddressing() const
{
    faceWeightsPtr_ = new FieldField<Field, scalar>(toPatch_.size());
    FieldField<Field, scalar>& faceWeights = *faceWeightsPtr_;

    faceDistancePtr_ = new scalarField(toPatch_.size(), GREAT);
    scalarField& faceDistance = *faceDistancePtr_;

    if (debug)
    {
        Info<< "projecting face centres" << endl;
    }

    const pointField& fromPatchPoints = fromPatch_.points();
    const typename FromPatch::FaceListType& fromPatchFaces = fromPatch_;
    const labelListList& fromPatchFaceFaces = fromPatch_.faceFaces();

    vectorField fromPatchFaceCentres(fromPatchFaces.size());

    forAll(fromPatchFaceCentres, faceI)
    {
        fromPatchFaceCentres[faceI] =
            fromPatchFaces[faceI].centre(fromPatchPoints);
    }

    const pointField& toPatchPoints = toPatch_.points();
    const typename ToPatch::FaceListType& toPatchFaces = toPatch_;

    const vectorField& projectionDirection = toPatch_.faceNormals();

    List<objectHit> proj =
        toPatch_.projectFaceCentres
        (
            fromPatch_,
            projectionDirection,
            alg_,
            dir_
        );

    faceAddressingPtr_ = new labelList(proj.size(), -1);
    labelList& faceAddressing = *faceAddressingPtr_;

    forAll(faceAddressing, faceI)
    {
        if (proj[faceI].hit())
        {
            // A hit exists
            faceAddressing[faceI] = proj[faceI].hitObject();

            const typename FromPatch::FaceType& hitFace =
                fromPatchFaces[faceAddressing[faceI]];

            pointHit curHit =
                hitFace.ray
                (
                    toPatchFaces[faceI].centre(toPatchPoints),
                    projectionDirection[faceI],
                    fromPatchPoints,
                    alg_,
                    dir_
                );

            // grab distance to target
            faceDistance[faceI] = curHit.distance();

            // grab face centre of the hit face
            const point& hitFaceCentre =
                fromPatchFaceCentres[faceAddressing[faceI]];

            // grab neighbours of hit face
            const labelList& neighbours =
                fromPatchFaceFaces[faceAddressing[faceI]];

            scalar m = mag(curHit.hitPoint() - hitFaceCentre);

            if
            (
                m < directHitTol                            // Direct hit
             || neighbours.empty()
            )
            {
                faceWeights.set(faceI, new scalarField(1));
                faceWeights[faceI][0] = 1.0;
            }
            else
            {
                // set interpolation faceWeights

                // The first coefficient corresponds to the centre face.
                // The rest is ordered in the same way as the faceFaces list.
                faceWeights.set(faceI, new scalarField(neighbours.size() + 1));

                faceWeights[faceI][0] = 1.0/m;

                forAll(neighbours, nI)
                {
                    faceWeights[faceI][nI + 1] =
                    1.0/
                    (
                        mag
                        (
                            fromPatchFaceCentres[neighbours[nI]]
                          - curHit.hitPoint()
                        )
                      + VSMALL
                    );
                }
            }

            faceWeights[faceI] /= sum(faceWeights[faceI]);
        }
        else
        {
            faceWeights.set(faceI, new scalarField(0));
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
