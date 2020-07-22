/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "PairCollision.H"
#include "PairModel.H"
#include "WallModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::PairCollision<CloudType>::cosPhiMinFlatWall = 1 - small;

template<class CloudType>
Foam::scalar Foam::PairCollision<CloudType>::flatWallDuplicateExclusion =
    sqrt(3*small);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::PairCollision<CloudType>::preInteraction()
{
    // Set accumulated quantities to zero
    forAllIter(typename CloudType, this->owner(), iter)
    {
        typename CloudType::parcelType& p = iter();

        p.f() = Zero;

        p.torque() = Zero;
    }
}


template<class CloudType>
void Foam::PairCollision<CloudType>::parcelInteraction()
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    label startOfRequests = Pstream::nRequests();

    il_.sendReferredData(this->owner().cellOccupancy(), pBufs);

    realRealInteraction();

    il_.receiveReferredData(pBufs, startOfRequests);

    realReferredInteraction();
}


template<class CloudType>
void Foam::PairCollision<CloudType>::realRealInteraction()
{
    // Direct interaction list (dil)
    const labelListList& dil = il_.dil();

    typename CloudType::parcelType* pA_ptr = nullptr;
    typename CloudType::parcelType* pB_ptr = nullptr;

    List<DynamicList<typename CloudType::parcelType*>>& cellOccupancy =
        this->owner().cellOccupancy();

    forAll(dil, realCelli)
    {
        // Loop over all Parcels in cell A (a)
        forAll(cellOccupancy[realCelli], a)
        {
            pA_ptr = cellOccupancy[realCelli][a];

            forAll(dil[realCelli], interactingCells)
            {
                List<typename CloudType::parcelType*> cellBParcels =
                    cellOccupancy[dil[realCelli][interactingCells]];

                // Loop over all Parcels in cell B (b)
                forAll(cellBParcels, b)
                {
                    pB_ptr = cellBParcels[b];

                    evaluatePair(*pA_ptr, *pB_ptr);
                }
            }

            // Loop over the other Parcels in cell A (aO)
            forAll(cellOccupancy[realCelli], aO)
            {
                pB_ptr = cellOccupancy[realCelli][aO];

                // Do not double-evaluate, compare pointers, arbitrary
                // order
                if (pB_ptr > pA_ptr)
                {
                    evaluatePair(*pA_ptr, *pB_ptr);
                }
            }
        }
    }
}


template<class CloudType>
void Foam::PairCollision<CloudType>::realReferredInteraction()
{
    // Referred interaction list (ril)
    const labelListList& ril = il_.ril();

    List<IDLList<typename CloudType::parcelType>>& referredParticles =
        il_.referredParticles();

    List<DynamicList<typename CloudType::parcelType*>>& cellOccupancy =
        this->owner().cellOccupancy();

    // Loop over all referred cells
    forAll(ril, refCelli)
    {
        IDLList<typename CloudType::parcelType>& refCellRefParticles =
            referredParticles[refCelli];

        const labelList& realCells = ril[refCelli];

        // Loop over all referred parcels in the referred cell

        forAllIter
        (
            typename IDLList<typename CloudType::parcelType>,
            refCellRefParticles,
            referredParcel
        )
        {
            // Loop over all real cells in that the referred cell is
            // to supply interactions to

            forAll(realCells, realCelli)
            {
                List<typename CloudType::parcelType*> realCellParcels =
                    cellOccupancy[realCells[realCelli]];

                forAll(realCellParcels, realParcelI)
                {
                    evaluatePair
                    (
                        *realCellParcels[realParcelI],
                        referredParcel()
                    );
                }
            }
        }
    }
}


template<class CloudType>
void Foam::PairCollision<CloudType>::wallInteraction()
{
    const polyMesh& mesh = this->owner().mesh();

    const labelListList& dil = il_.dil();

    const labelListList& directWallFaces = il_.dwfil();

    const labelList& patchID = mesh.boundaryMesh().patchID();

    const volVectorField& U = mesh.lookupObject<volVectorField>(il_.UName());

    List<DynamicList<typename CloudType::parcelType*>>& cellOccupancy =
        this->owner().cellOccupancy();

    // Storage for the wall interaction sites
    DynamicList<point> flatSitePoints;
    DynamicList<scalar> flatSiteExclusionDistancesSqr;
    DynamicList<WallSiteData<vector>> flatSiteData;
    DynamicList<point> otherSitePoints;
    DynamicList<scalar> otherSiteDistances;
    DynamicList<WallSiteData<vector>> otherSiteData;
    DynamicList<point> sharpSitePoints;
    DynamicList<scalar> sharpSiteExclusionDistancesSqr;
    DynamicList<WallSiteData<vector>> sharpSiteData;

    forAll(dil, realCelli)
    {
        // The real wall faces in range of this real cell
        const labelList& realWallFaces = directWallFaces[realCelli];

        // Loop over all Parcels in cell
        forAll(cellOccupancy[realCelli], cellParticleI)
        {
            flatSitePoints.clear();
            flatSiteExclusionDistancesSqr.clear();
            flatSiteData.clear();
            otherSitePoints.clear();
            otherSiteDistances.clear();
            otherSiteData.clear();
            sharpSitePoints.clear();
            sharpSiteExclusionDistancesSqr.clear();
            sharpSiteData.clear();

            typename CloudType::parcelType& p =
                *cellOccupancy[realCelli][cellParticleI];

            const point& pos = p.position();

            scalar r = wallModel_->pREff(p);

            // real wallFace interactions

            forAll(realWallFaces, realWallFacei)
            {
                label realFacei = realWallFaces[realWallFacei];

                pointHit nearest = mesh.faces()[realFacei].nearestPoint
                (
                    pos,
                    mesh.points()
                );

                if (nearest.distance() < r)
                {
                    vector normal = mesh.faceAreas()[realFacei];

                    normal /= mag(normal);

                    const vector& nearPt = nearest.rawPoint();

                    vector pW = nearPt - pos;

                    scalar normalAlignment = normal & pW/(mag(pW) + small);

                    // Find the patchIndex and wallData for WallSiteData object
                    label patchi = patchID[realFacei - mesh.nInternalFaces()];

                    label patchFacei =
                        realFacei - mesh.boundaryMesh()[patchi].start();

                    WallSiteData<vector> wSD
                    (
                        patchi,
                        U.boundaryField()[patchi][patchFacei]
                    );

                    if (normalAlignment > cosPhiMinFlatWall)
                    {
                        // Guard against a flat interaction being
                        // present on the boundary of two or more
                        // faces, which would create duplicate contact
                        // points. Duplicates are discarded.
                        if
                        (
                            !duplicatePointInList
                            (
                                flatSitePoints,
                                nearPt,
                                sqr(r*flatWallDuplicateExclusion)
                            )
                        )
                        {
                            flatSitePoints.append(nearPt);

                            flatSiteExclusionDistancesSqr.append
                            (
                                sqr(r) - sqr(nearest.distance())
                            );

                            flatSiteData.append(wSD);
                        }
                    }
                    else
                    {
                        otherSitePoints.append(nearPt);

                        otherSiteDistances.append(nearest.distance());

                        otherSiteData.append(wSD);
                    }
                }
            }

            // referred wallFace interactions

            // The labels of referred wall faces in range of this real cell
            const labelList& cellRefWallFaces = il_.rwfilInverse()[realCelli];

            forAll(cellRefWallFaces, rWFI)
            {
                label refWallFacei = cellRefWallFaces[rWFI];

                const referredWallFace& rwf =
                    il_.referredWallFaces()[refWallFacei];

                const pointField& pts = rwf.points();

                pointHit nearest = rwf.nearestPoint(pos, pts);

                if (nearest.distance() < r)
                {
                    const vector normal = rwf.normal(pts);
                    const vector& nearPt = nearest.rawPoint();

                    vector pW = nearPt - pos;

                    scalar normalAlignment = normal & pW/mag(pW);

                    // Find the patchIndex and wallData for WallSiteData object

                    WallSiteData<vector> wSD
                    (
                        rwf.patchIndex(),
                        il_.referredWallData()[refWallFacei]
                    );

                    if (normalAlignment > cosPhiMinFlatWall)
                    {
                        // Guard against a flat interaction being
                        // present on the boundary of two or more
                        // faces, which would create duplicate contact
                        // points. Duplicates are discarded.
                        if
                        (
                            !duplicatePointInList
                            (
                                flatSitePoints,
                                nearPt,
                                sqr(r*flatWallDuplicateExclusion)
                            )
                        )
                        {
                            flatSitePoints.append(nearPt);

                            flatSiteExclusionDistancesSqr.append
                            (
                                sqr(r) - sqr(nearest.distance())
                            );

                            flatSiteData.append(wSD);
                        }
                    }
                    else
                    {
                        otherSitePoints.append(nearPt);

                        otherSiteDistances.append(nearest.distance());

                        otherSiteData.append(wSD);
                    }
                }
            }

            // All flat interaction sites found, now classify the
            // other sites as being in range of a flat interaction, or
            // a sharp interaction, being aware of not duplicating the
            // sharp interaction sites.

            // The "other" sites need to evaluated in order of
            // ascending distance to their nearest point so that
            // grouping occurs around the closest in any group

            labelList sortedOtherSiteIndices;

            sortedOrder(otherSiteDistances, sortedOtherSiteIndices);

            forAll(sortedOtherSiteIndices, siteI)
            {
                label orderedIndex = sortedOtherSiteIndices[siteI];

                const point& otherPt = otherSitePoints[orderedIndex];

                if
                (
                    !duplicatePointInList
                    (
                        flatSitePoints,
                        otherPt,
                        flatSiteExclusionDistancesSqr
                    )
                )
                {
                    // Not in range of a flat interaction, must be a
                    // sharp interaction.

                    if
                    (
                        !duplicatePointInList
                        (
                            sharpSitePoints,
                            otherPt,
                            sharpSiteExclusionDistancesSqr
                        )
                    )
                    {
                        sharpSitePoints.append(otherPt);

                        sharpSiteExclusionDistancesSqr.append
                        (
                            sqr(r) - sqr(otherSiteDistances[orderedIndex])
                        );

                        sharpSiteData.append(otherSiteData[orderedIndex]);
                    }
                }
            }

            evaluateWall
            (
                p,
                flatSitePoints,
                flatSiteData,
                sharpSitePoints,
                sharpSiteData
            );
        }
    }
}


template<class CloudType>
bool Foam::PairCollision<CloudType>::duplicatePointInList
(
    const DynamicList<point>& existingPoints,
    const point& pointToTest,
    scalar duplicateRangeSqr
) const
{
    forAll(existingPoints, i)
    {
        if (magSqr(existingPoints[i] - pointToTest) < duplicateRangeSqr)
        {
            return true;
        }
    }

    return false;
}


template<class CloudType>
bool Foam::PairCollision<CloudType>::duplicatePointInList
(
    const DynamicList<point>& existingPoints,
    const point& pointToTest,
    const scalarList& duplicateRangeSqr
) const
{
    forAll(existingPoints, i)
    {
        if (magSqr(existingPoints[i] - pointToTest) < duplicateRangeSqr[i])
        {
            return true;
        }
    }

    return false;
}


template<class CloudType>
void Foam::PairCollision<CloudType>::postInteraction()
{
    // Delete any collision records where no collision occurred this step

    forAllIter(typename CloudType, this->owner(), iter)
    {
        typename CloudType::parcelType& p = iter();

        p.collisionRecords().update();
    }
}


template<class CloudType>
void Foam::PairCollision<CloudType>::evaluatePair
(
    typename CloudType::parcelType& pA,
    typename CloudType::parcelType& pB
) const
{
    pairModel_->evaluatePair(pA, pB);
}


template<class CloudType>
void Foam::PairCollision<CloudType>::evaluateWall
(
    typename CloudType::parcelType& p,
    const List<point>& flatSitePoints,
    const List<WallSiteData<vector>>& flatSiteData,
    const List<point>& sharpSitePoints,
    const List<WallSiteData<vector>>& sharpSiteData
) const
{
    wallModel_->evaluateWall
    (
        p,
        flatSitePoints,
        flatSiteData,
        sharpSitePoints,
        sharpSiteData
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PairCollision<CloudType>::PairCollision
(
    const dictionary& dict,
    CloudType& owner
)
:
    CollisionModel<CloudType>(dict, owner, typeName),
    pairModel_
    (
        PairModel<CloudType>::New
        (
            this->coeffDict(),
            this->owner()
        )
    ),
    wallModel_
    (
        WallModel<CloudType>::New
        (
            this->coeffDict(),
            this->owner()
        )
    ),
    il_
    (
        owner.mesh(),
        this->coeffDict().template lookup<scalar>("maxInteractionDistance"),
        Switch
        (
            this->coeffDict().lookupOrDefault
            (
                "writeReferredParticleCloud",
                false
            )
        ),
        this->coeffDict().lookupOrDefault("U", word("U"))
    )
{}


template<class CloudType>
Foam::PairCollision<CloudType>::PairCollision
(
    const PairCollision<CloudType>& cm
)
:
    CollisionModel<CloudType>(cm),
    pairModel_(nullptr),
    wallModel_(nullptr),
    il_(cm.owner().mesh())
{
    // Need to clone to PairModel and WallModel
    NotImplemented;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PairCollision<CloudType>::~PairCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::PairCollision<CloudType>::nSubCycles() const
{
    label nSubCycles = 1;

    if (pairModel_->controlsTimestep())
    {
        label nPairSubCycles = returnReduce
        (
            pairModel_->nSubCycles(), maxOp<label>()
        );

        nSubCycles = max(nSubCycles, nPairSubCycles);
    }

    if (wallModel_->controlsTimestep())
    {
        label nWallSubCycles = returnReduce
        (
            wallModel_->nSubCycles(), maxOp<label>()
        );

        nSubCycles = max(nSubCycles, nWallSubCycles);
    }

    return nSubCycles;
}


template<class CloudType>
void Foam::PairCollision<CloudType>::collide()
{
    preInteraction();

    parcelInteraction();

    wallInteraction();

    postInteraction();
}


// ************************************************************************* //
