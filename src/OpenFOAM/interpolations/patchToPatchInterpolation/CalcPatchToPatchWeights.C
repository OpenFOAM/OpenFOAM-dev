/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

    pointDistancePtr_ = new scalarField(toPatch_.nPoints(), great);
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

    forAll(pointAddressing, pointi)
    {
        doWeights = false;

        const typename FromPatch::FaceType& hitFace =
            fromPatchFaces[proj[pointi].hitObject()];

        point hitPoint = Zero;

        if (proj[pointi].hit())
        {
            // A hit exists
            doWeights = true;

            pointAddressing[pointi] = proj[pointi].hitObject();

            pointHit curHit =
                hitFace.ray
                (
                    toPatchPoints[pointi],
                    projectionDirection[pointi],
                    fromPatchPoints,
                    alg_,
                    dir_
                );

            // Grab distance to target
            if (dir_ == intersection::CONTACT_SPHERE)
            {
                pointDistance[pointi] =
                    hitFace.contactSphereDiameter
                    (
                        toPatchPoints[pointi],
                        projectionDirection[pointi],
                        fromPatchPoints
                    );
            }
            else
            {
                pointDistance[pointi] = curHit.distance();
            }

            // Grab hit point
            hitPoint = curHit.hitPoint();
        }
        else if (projectionTol_ > small)
        {
            // Check for a near miss
            pointHit ph =
                hitFace.ray
                (
                    toPatchPoints[pointi],
                    projectionDirection[pointi],
                    fromPatchPoints,
                    alg_,
                    dir_
                );

            scalar dist =
                Foam::mag
                (
                    toPatchPoints[pointi]
                  + projectionDirection[pointi]*ph.distance()
                  - ph.missPoint()
                );

            // Calculate the local tolerance
            scalar minEdgeLength = great;

            // Do shortest edge of hit object
            edgeList hitFaceEdges =
                fromPatchFaces[proj[pointi].hitObject()].edges();

            forAll(hitFaceEdges, edgeI)
            {
                minEdgeLength =
                    min
                    (
                        minEdgeLength,
                        hitFaceEdges[edgeI].mag(fromPatchPoints)
                    );
            }

            const labelList& curEdges = toPatchPointEdges[pointi];

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

                pointAddressing[pointi] = proj[pointi].hitObject();

                // Grab nearest point on face as hit point
                hitPoint = ph.missPoint();

                // Grab distance to target
                if (dir_ == intersection::CONTACT_SPHERE)
                {
                    pointDistance[pointi] =
                        hitFace.contactSphereDiameter
                        (
                            toPatchPoints[pointi],
                            projectionDirection[pointi],
                            fromPatchPoints
                        );
                }
                else
                {
                    pointDistance[pointi] =
                        (
                            projectionDirection[pointi]
                            /mag(projectionDirection[pointi])
                        )
                      & (hitPoint - toPatchPoints[pointi]);
                }
            }
        }

        if (doWeights)
        {
            // Set interpolation pointWeights
            pointWeights.set(pointi, new scalarField(hitFace.size()));

            pointField hitFacePoints = hitFace.points(fromPatchPoints);

            forAll(hitFacePoints, masterPointi)
            {
                pointWeights[pointi][masterPointi] =
                    1.0/
                    (
                        mag
                        (
                            hitFacePoints[masterPointi]
                          - hitPoint
                        )
                      + vSmall
                    );
            }

            pointWeights[pointi] /= sum(pointWeights[pointi]);
        }
        else
        {
            pointWeights.set(pointi, new scalarField(0));
        }
    }
}


template<class FromPatch, class ToPatch>
void PatchToPatchInterpolation<FromPatch, ToPatch>::calcFaceAddressing() const
{
    faceWeightsPtr_ = new FieldField<Field, scalar>(toPatch_.size());
    FieldField<Field, scalar>& faceWeights = *faceWeightsPtr_;

    faceDistancePtr_ = new scalarField(toPatch_.size(), great);
    scalarField& faceDistance = *faceDistancePtr_;

    if (debug)
    {
        Info<< "projecting face centres" << endl;
    }

    const pointField& fromPatchPoints = fromPatch_.points();
    const typename FromPatch::FaceListType& fromPatchFaces = fromPatch_;
    const labelListList& fromPatchFaceFaces = fromPatch_.faceFaces();

    vectorField fromPatchFaceCentres(fromPatchFaces.size());

    forAll(fromPatchFaceCentres, facei)
    {
        fromPatchFaceCentres[facei] =
            fromPatchFaces[facei].centre(fromPatchPoints);
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

    forAll(faceAddressing, facei)
    {
        if (proj[facei].hit())
        {
            // A hit exists
            faceAddressing[facei] = proj[facei].hitObject();

            const typename FromPatch::FaceType& hitFace =
                fromPatchFaces[faceAddressing[facei]];

            pointHit curHit =
                hitFace.ray
                (
                    toPatchFaces[facei].centre(toPatchPoints),
                    projectionDirection[facei],
                    fromPatchPoints,
                    alg_,
                    dir_
                );

            // grab distance to target
            faceDistance[facei] = curHit.distance();

            // grab face centre of the hit face
            const point& hitFaceCentre =
                fromPatchFaceCentres[faceAddressing[facei]];

            // grab neighbours of hit face
            const labelList& neighbours =
                fromPatchFaceFaces[faceAddressing[facei]];

            scalar m = mag(curHit.hitPoint() - hitFaceCentre);

            if
            (
                m < directHitTol                            // Direct hit
             || neighbours.empty()
            )
            {
                faceWeights.set(facei, new scalarField(1));
                faceWeights[facei][0] = 1.0;
            }
            else
            {
                // set interpolation faceWeights

                // The first coefficient corresponds to the centre face.
                // The rest is ordered in the same way as the faceFaces list.
                faceWeights.set(facei, new scalarField(neighbours.size() + 1));

                faceWeights[facei][0] = 1.0/m;

                forAll(neighbours, nI)
                {
                    faceWeights[facei][nI + 1] =
                    1.0/
                    (
                        mag
                        (
                            fromPatchFaceCentres[neighbours[nI]]
                          - curHit.hitPoint()
                        )
                      + vSmall
                    );
                }
            }

            faceWeights[facei] /= sum(faceWeights[facei]);
        }
        else
        {
            faceWeights.set(facei, new scalarField(0));
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
