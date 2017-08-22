/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "IOPosition.H"

#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "cyclicACMIPolyPatch.H"
#include "processorPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wallPolyPatch.H"
#include "wedgePolyPatch.H"
#include "meshTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class TrackCloudType>
void Foam::particle::prepareForParallelTransfer
(
    const label patchi,
    TrackCloudType& cloud,
    trackingData& td
)
{
    // Convert the face index to be local to the processor patch
    facei_ = mesh_.boundaryMesh()[patchi].whichFace(facei_);
}


template<class TrackCloudType>
void Foam::particle::correctAfterParallelTransfer
(
    const label patchi,
    TrackCloudType& cloud,
    trackingData& td
)
{
    const coupledPolyPatch& ppp =
        refCast<const coupledPolyPatch>(mesh_.boundaryMesh()[patchi]);

    if (!ppp.parallel())
    {
        const tensor& T =
        (
            ppp.forwardT().size() == 1
          ? ppp.forwardT()[0]
          : ppp.forwardT()[facei_]
        );
        transformProperties(T);
    }
    else if (ppp.separated())
    {
        const vector& s =
        (
            (ppp.separation().size() == 1)
          ? ppp.separation()[0]
          : ppp.separation()[facei_]
        );
        transformProperties(-s);
    }

    // Set the topology
    celli_ = ppp.faceCells()[facei_];
    facei_ += ppp.start();
    tetFacei_ = facei_;
    // Faces either side of a coupled patch are numbered in opposite directions
    // as their normals both point away from their connected cells. The tet
    // point therefore counts in the opposite direction from the base point.
    tetPti_ = mesh_.faces()[tetFacei_].size() - 1 - tetPti_;

    // Reflect to account for the change of triangle orientation in the new cell
    reflect();

    // Note that the position does not need transforming explicitly. The face-
    // triangle on the receive patch is the transformation of the one on the
    // send patch, so whilst the barycentric coordinates remain the same, the
    // change of triangle implicitly transforms the position.
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackCloudType>
void Foam::particle::readFields(TrackCloudType& c)
{
    bool valid = c.size();

    IOobject procIO(c.fieldIOobject("origProcId", IOobject::MUST_READ));

    bool haveFile = procIO.typeHeaderOk<IOField<label>>(true);

    IOField<label> origProcId(procIO, valid && haveFile);
    c.checkFieldIOobject(c, origProcId);
    IOField<label> origId
    (
        c.fieldIOobject("origId", IOobject::MUST_READ),
        valid && haveFile
    );
    c.checkFieldIOobject(c, origId);

    label i = 0;
    forAllIter(typename TrackCloudType, c, iter)
    {
        particle& p = iter();

        p.origProc_ = origProcId[i];
        p.origId_ = origId[i];
        i++;
    }
}


template<class TrackCloudType>
void Foam::particle::writeFields(const TrackCloudType& c)
{
    label np = c.size();

    IOPosition<TrackCloudType> ioP(c);
    ioP.write(np > 0);

    IOField<label> origProc
    (
        c.fieldIOobject("origProcId", IOobject::NO_READ),
        np
    );
    IOField<label> origId
    (
        c.fieldIOobject("origId", IOobject::NO_READ),
        np
    );

    label i = 0;
    forAllConstIter(typename TrackCloudType, c, iter)
    {
        origProc[i] = iter().origProc_;
        origId[i] = iter().origId_;
        i++;
    }

    origProc.write(np > 0);
    origId.write(np > 0);
}


template<class TrackCloudType>
void Foam::particle::hitFace
(
    const vector& direction,
    TrackCloudType& cloud,
    trackingData& td
)
{
    typename TrackCloudType::particleType& p =
        static_cast<typename TrackCloudType::particleType&>(*this);
    typename TrackCloudType::particleType::trackingData& ttd =
        static_cast<typename TrackCloudType::particleType::trackingData&>(td);

    if (!onFace())
    {
        return;
    }
    else if (onInternalFace())
    {
        changeCell();
    }
    else if (onBoundaryFace())
    {
        const tetIndices faceHitTetIs(celli_, tetFacei_, tetPti_);

        if
        (
           !p.hitPatch
            (
                mesh_.boundaryMesh()[patch()],
                cloud,
                ttd,
                patch(),
                stepFraction(),
                faceHitTetIs
            )
        )
        {
            const polyPatch& patch = mesh_.boundaryMesh()[this->patch()];

            if (isA<wedgePolyPatch>(patch))
            {
                p.hitWedgePatch
                (
                    static_cast<const wedgePolyPatch&>(patch), cloud, ttd
                );
            }
            else if (isA<symmetryPlanePolyPatch>(patch))
            {
                p.hitSymmetryPlanePatch
                (
                    static_cast<const symmetryPlanePolyPatch&>(patch),
                    cloud,
                    ttd
                );
            }
            else if (isA<symmetryPolyPatch>(patch))
            {
                p.hitSymmetryPatch
                (
                    static_cast<const symmetryPolyPatch&>(patch), cloud, ttd
                );
            }
            else if (isA<cyclicPolyPatch>(patch))
            {
                p.hitCyclicPatch
                (
                    static_cast<const cyclicPolyPatch&>(patch), cloud, ttd
                );
            }
            else if (isA<cyclicACMIPolyPatch>(patch))
            {
                p.hitCyclicACMIPatch
                (
                    static_cast<const cyclicACMIPolyPatch&>(patch),
                    cloud,
                    ttd,
                    direction
                );
            }
            else if (isA<cyclicAMIPolyPatch>(patch))
            {
                p.hitCyclicAMIPatch
                (
                    static_cast<const cyclicAMIPolyPatch&>(patch),
                    cloud,
                    ttd,
                    direction
                );
            }
            else if (isA<processorPolyPatch>(patch))
            {
                p.hitProcessorPatch
                (
                    static_cast<const processorPolyPatch&>(patch), cloud, ttd
                );
            }
            else if (isA<wallPolyPatch>(patch))
            {
                p.hitWallPatch
                (
                    static_cast<const wallPolyPatch&>(patch),
                    cloud,
                    ttd,
                    faceHitTetIs
                );
            }
            else
            {
                p.hitPatch(patch, cloud, ttd);
            }
        }
    }
}


template<class TrackCloudType>
void Foam::particle::trackToAndHitFace
(
    const vector& direction,
    const scalar fraction,
    TrackCloudType& cloud,
    trackingData& td
)
{
    trackToFace(direction, fraction);

    if (onBoundaryFace())
    {
        changeToMasterPatch();
    }

    hitFace(direction, cloud, td);
}


template<class TrackCloudType>
bool Foam::particle::hitPatch
(
    const polyPatch&,
    TrackCloudType&,
    trackingData&,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}


template<class TrackCloudType>
void Foam::particle::hitWedgePatch
(
    const wedgePolyPatch& wpp,
    TrackCloudType&,
    trackingData&
)
{
    FatalErrorInFunction
        << "Hitting a wedge patch should not be possible."
        << abort(FatalError);

    vector nf = normal();
    nf /= mag(nf);

    transformProperties(I - 2.0*nf*nf);
}


template<class TrackCloudType>
void Foam::particle::hitSymmetryPlanePatch
(
    const symmetryPlanePolyPatch& spp,
    TrackCloudType&,
    trackingData&
)
{
    vector nf = normal();
    nf /= mag(nf);

    transformProperties(I - 2.0*nf*nf);
}


template<class TrackCloudType>
void Foam::particle::hitSymmetryPatch
(
    const symmetryPolyPatch& spp,
    TrackCloudType&,
    trackingData&
)
{
    vector nf = normal();
    nf /= mag(nf);

    transformProperties(I - 2.0*nf*nf);
}


template<class TrackCloudType>
void Foam::particle::hitCyclicPatch
(
    const cyclicPolyPatch& cpp,
    TrackCloudType&,
    trackingData&
)
{
    const cyclicPolyPatch& receiveCpp = cpp.neighbPatch();
    const label receiveFacei = receiveCpp.whichFace(facei_);

    // Set the topology
    facei_ = tetFacei_ = cpp.transformGlobalFace(facei_);
    celli_ = mesh_.faceOwner()[facei_];
    // See note in correctAfterParallelTransfer for tetPti addressing ...
    tetPti_ = mesh_.faces()[tetFacei_].size() - 1 - tetPti_;

    // Reflect to account for the change of triangle orientation in the new cell
    reflect();

    // Transform the properties
    if (!receiveCpp.parallel())
    {
        const tensor& T =
        (
            receiveCpp.forwardT().size() == 1
          ? receiveCpp.forwardT()[0]
          : receiveCpp.forwardT()[receiveFacei]
        );
        transformProperties(T);
    }
    else if (receiveCpp.separated())
    {
        const vector& s =
        (
            (receiveCpp.separation().size() == 1)
          ? receiveCpp.separation()[0]
          : receiveCpp.separation()[receiveFacei]
        );
        transformProperties(-s);
    }
}


template<class TrackCloudType>
void Foam::particle::hitCyclicAMIPatch
(
    const cyclicAMIPolyPatch& cpp,
    TrackCloudType& cloud,
    trackingData& td,
    const vector& direction
)
{
    vector pos = position();

    const cyclicAMIPolyPatch& receiveCpp = cpp.neighbPatch();
    const label sendFacei = cpp.whichFace(facei_);
    const label receiveFacei = cpp.pointFace(sendFacei, direction, pos);

    if (receiveFacei < 0)
    {
        // If the patch face of the particle is not known assume that the
        // particle is lost and mark it to be deleted.
        td.keepParticle = false;
        WarningInFunction
            << "Particle lost across " << cyclicAMIPolyPatch::typeName
            << " patches " << cpp.name() << " and " << receiveCpp.name()
            << " at position " << pos << endl;
    }

    // Set the topology
    facei_ = tetFacei_ = receiveFacei + receiveCpp.start();

    // Locate the particle on the recieving side
    vector directionT = direction;
    cpp.reverseTransformDirection(directionT, sendFacei);
    locate
    (
        pos,
        &directionT,
        mesh_.faceOwner()[facei_],
        false,
        "Particle crossed between " + cyclicAMIPolyPatch::typeName +
        " patches " + cpp.name() + " and " + receiveCpp.name() +
        " to a location outside of the mesh."
    );

    // The particle must remain associated with a face for the tracking to
    // register as incomplete
    facei_ = tetFacei_;

    // Transform the properties
    if (!receiveCpp.parallel())
    {
        const tensor& T =
        (
            receiveCpp.forwardT().size() == 1
          ? receiveCpp.forwardT()[0]
          : receiveCpp.forwardT()[receiveFacei]
        );
        transformProperties(T);
    }
    else if (receiveCpp.separated())
    {
        const vector& s =
        (
            (receiveCpp.separation().size() == 1)
          ? receiveCpp.separation()[0]
          : receiveCpp.separation()[receiveFacei]
        );
        transformProperties(-s);
    }
}


template<class TrackCloudType>
void Foam::particle::hitCyclicACMIPatch
(
    const cyclicACMIPolyPatch& cpp,
    TrackCloudType& cloud,
    trackingData& td,
    const vector& direction
)
{
    const label localFacei = cpp.whichFace(facei_);

    // If the mask is within the patch tolerance at either end, then we can
    // assume an interaction with the appropriate part of the ACMI pair.
    const scalar mask = cpp.mask()[localFacei];
    bool couple = mask >= 1 - cpp.tolerance();
    bool nonOverlap = mask <= cpp.tolerance();

    // If the mask is an intermediate value, then we search for a location on
    // the other side of the AMI. If we can't find a location, then we assume
    // that we have hit the non-overlap patch.
    if (!couple && !nonOverlap)
    {
        vector pos = position();
        couple = cpp.pointFace(localFacei, direction, pos) >= 0;
        nonOverlap = !couple;
    }

    if (couple)
    {
        hitCyclicAMIPatch(cpp, cloud, td, direction);
    }
    else
    {
        // Move to the face associated with the non-overlap patch and redo the
        // face interaction.
        tetFacei_ = facei_ = cpp.nonOverlapPatch().start() + localFacei;
        hitFace(direction, cloud, td);
    }
}


template<class TrackCloudType>
void Foam::particle::hitProcessorPatch
(
    const processorPolyPatch&,
    TrackCloudType&,
    trackingData&
)
{}


template<class TrackCloudType>
void Foam::particle::hitWallPatch
(
    const wallPolyPatch&,
    TrackCloudType&,
    trackingData&,
    const tetIndices&
)
{}


template<class TrackCloudType>
void Foam::particle::hitPatch(const polyPatch&, TrackCloudType&, trackingData&)
{}


// ************************************************************************* //
