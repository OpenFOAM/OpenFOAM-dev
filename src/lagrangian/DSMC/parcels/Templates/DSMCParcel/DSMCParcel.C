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

#include "DSMCParcel.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
bool Foam::DSMCParcel<ParcelType>::move
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = cloud.pMesh();

    // For reduced-D cases, the velocity used to track needs to be
    // constrained, but the actual U_ of the parcel must not be
    // altered or used, as it is altered by patch interactions an
    // needs to retain its 3D value for collision purposes.
    vector Utracking = U_;

    while (td.keepParticle && !td.switchProcessor && p.stepFraction() < 1)
    {
        Utracking = U_;

        // Apply correction to velocity to constrain tracking for
        // reduced-D cases
        meshTools::constrainDirection(mesh, mesh.solutionD(), Utracking);

        // Deviation from the mesh centre for reduced-D cases
        const vector d = p.deviationFromMeshCentre();

        const scalar f = 1 - p.stepFraction();
        p.trackToAndHitFace(f*trackTime*Utracking - d, f, cloud, td);
    }

    return td.keepParticle;
}


template<class ParcelType>
template<class TrackCloudType>
bool Foam::DSMCParcel<ParcelType>::hitPatch(TrackCloudType&, trackingData&)
{
    return false;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::DSMCParcel<ParcelType>::hitProcessorPatch
(
    TrackCloudType&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::DSMCParcel<ParcelType>::hitWallPatch
(
    TrackCloudType& cloud,
    trackingData&
)
{
    const label wppIndex = this->patch();

    const wallPolyPatch& wpp =
        static_cast<const wallPolyPatch&>
        (
            this->mesh().boundaryMesh()[wppIndex]
        );

    const label wppLocalFace = wpp.whichFace(this->face());

    const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);

    const scalar deltaT = cloud.pMesh().time().deltaTValue();

    const constantProperties& constProps(cloud.constProps(typeId_));

    scalar m = constProps.mass();

    vector nw = wpp.faceAreas()[wppLocalFace];
    nw /= mag(nw);

    scalar U_dot_nw = U_ & nw;

    vector Ut = U_ - U_dot_nw*nw;

    scalar invMagUnfA = 1/max(mag(U_dot_nw)*fA, vSmall);

    cloud.rhoNBF()[wppIndex][wppLocalFace] += invMagUnfA;

    cloud.rhoMBF()[wppIndex][wppLocalFace] += m*invMagUnfA;

    cloud.linearKEBF()[wppIndex][wppLocalFace] +=
        0.5*m*(U_ & U_)*invMagUnfA;

    cloud.internalEBF()[wppIndex][wppLocalFace] += Ei_*invMagUnfA;

    cloud.iDofBF()[wppIndex][wppLocalFace] +=
        constProps.internalDegreesOfFreedom()*invMagUnfA;

    cloud.momentumBF()[wppIndex][wppLocalFace] += m*Ut*invMagUnfA;

    // pre-interaction energy
    scalar preIE = 0.5*m*(U_ & U_) + Ei_;

    // pre-interaction momentum
    vector preIMom = m*U_;

    cloud.wallInteraction().correct(*this);

    U_dot_nw = U_ & nw;

    Ut = U_ - U_dot_nw*nw;

    invMagUnfA = 1/max(mag(U_dot_nw)*fA, vSmall);

    cloud.rhoNBF()[wppIndex][wppLocalFace] += invMagUnfA;

    cloud.rhoMBF()[wppIndex][wppLocalFace] += m*invMagUnfA;

    cloud.linearKEBF()[wppIndex][wppLocalFace] +=
        0.5*m*(U_ & U_)*invMagUnfA;

    cloud.internalEBF()[wppIndex][wppLocalFace] += Ei_*invMagUnfA;

    cloud.iDofBF()[wppIndex][wppLocalFace] +=
        constProps.internalDegreesOfFreedom()*invMagUnfA;

    cloud.momentumBF()[wppIndex][wppLocalFace] += m*Ut*invMagUnfA;

    // post-interaction energy
    scalar postIE = 0.5*m*(U_ & U_) + Ei_;

    // post-interaction momentum
    vector postIMom = m*U_;

    scalar deltaQ = cloud.nParticle()*(preIE - postIE)/(deltaT*fA);

    vector deltaFD = cloud.nParticle()*(preIMom - postIMom)/(deltaT*fA);

    cloud.qBF()[wppIndex][wppLocalFace] += deltaQ;

    cloud.fDBF()[wppIndex][wppLocalFace] += deltaFD;
}


template<class ParcelType>
void Foam::DSMCParcel<ParcelType>::transformProperties(const tensor& T)
{
    ParcelType::transformProperties(T);
    U_ = transform(T, U_);
}


template<class ParcelType>
void Foam::DSMCParcel<ParcelType>::transformProperties
(
    const vector& separation
)
{
    ParcelType::transformProperties(separation);
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "DSMCParcelIO.C"

// ************************************************************************* //
