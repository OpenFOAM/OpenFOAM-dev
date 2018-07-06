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

#include "IOPosition.H"

#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "cyclicACMIPolyPatch.H"
#include "cyclicRepeatAMIPolyPatch.H"
#include "processorPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wallPolyPatch.H"
#include "wedgePolyPatch.H"
#include "meshTools.H"

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
        changeToMasterPatch();

        if (!p.hitPatch(cloud, ttd))
        {
            const polyPatch& patch = mesh_.boundaryMesh()[p.patch()];

            if (isA<wedgePolyPatch>(patch))
            {
                p.hitWedgePatch(cloud, ttd);
            }
            else if (isA<symmetryPlanePolyPatch>(patch))
            {
                p.hitSymmetryPlanePatch(cloud, ttd);
            }
            else if (isA<symmetryPolyPatch>(patch))
            {
                p.hitSymmetryPatch(cloud, ttd);
            }
            else if (isA<cyclicPolyPatch>(patch))
            {
                p.hitCyclicPatch(cloud, ttd);
            }
            else if (isA<cyclicACMIPolyPatch>(patch))
            {
                p.hitCyclicACMIPatch(cloud, ttd, direction);
            }
            else if (isA<cyclicAMIPolyPatch>(patch))
            {
                p.hitCyclicAMIPatch(cloud, ttd, direction);
            }
            else if (isA<cyclicRepeatAMIPolyPatch>(patch))
            {
                p.hitCyclicRepeatAMIPatch(cloud, ttd, direction);
            }
            else if (isA<processorPolyPatch>(patch))
            {
                p.hitProcessorPatch(cloud, ttd);
            }
            else if (isA<wallPolyPatch>(patch))
            {
                p.hitWallPatch(cloud, ttd);
            }
            else
            {
                td.keepParticle = false;
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

    hitFace(direction, cloud, td);
}


template<class TrackCloudType>
bool Foam::particle::hitPatch(TrackCloudType&, trackingData&)
{
    return false;
}


template<class TrackCloudType>
void Foam::particle::hitWedgePatch(TrackCloudType& cloud, trackingData& td)
{
    FatalErrorInFunction
        << "Hitting a wedge patch should not be possible."
        << abort(FatalError);

    hitSymmetryPatch(cloud, td);
}


template<class TrackCloudType>
void Foam::particle::hitSymmetryPlanePatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    hitSymmetryPatch(cloud, td);
}


template<class TrackCloudType>
void Foam::particle::hitSymmetryPatch(TrackCloudType&, trackingData&)
{
    const vector nf = normal();
    transformProperties(I - 2.0*nf*nf);
}


template<class TrackCloudType>
void Foam::particle::hitCyclicPatch(TrackCloudType&, trackingData&)
{
    const cyclicPolyPatch& cpp =
        static_cast<const cyclicPolyPatch&>(mesh_.boundaryMesh()[patch()]);
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
    TrackCloudType&,
    trackingData& td,
    const vector& direction
)
{
    vector pos = position();

    const cyclicAMIPolyPatch& cpp =
        static_cast<const cyclicAMIPolyPatch&>(mesh_.boundaryMesh()[patch()]);
    const cyclicAMIPolyPatch& receiveCpp = cpp.neighbPatch();
    const label sendFacei = cpp.whichFace(facei_);
    const labelPair receiveIs = cpp.pointAMIAndFace(sendFacei, direction, pos);
    const label receiveAMIi = receiveIs.first();
    const label receiveFacei = receiveIs.second();

    if (receiveFacei < 0)
    {
        // If the patch face of the particle is not known assume that the
        // particle is lost and mark it to be deleted.
        td.keepParticle = false;
        WarningInFunction
            << "Particle lost across " << cyclicAMIPolyPatch::typeName
            << " patches " << cpp.name() << " and " << receiveCpp.name()
            << " at position " << pos << endl;
        return;
    }

    // Set the topology
    facei_ = tetFacei_ = receiveFacei + receiveCpp.start();

    // Locate the particle on the receiving side
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
    const vectorTensorTransform& T =
        cpp.owner()
      ? cpp.AMITransforms()[receiveAMIi]
      : cpp.neighbPatch().AMITransforms()[receiveAMIi];
    if (T.hasR())
    {
        transformProperties(T.R());
    }
    else if (T.t() != vector::zero)
    {
        transformProperties(T.t());
    }
}


template<class TrackCloudType>
void Foam::particle::hitCyclicACMIPatch
(
    TrackCloudType& cloud,
    trackingData& td,
    const vector& direction
)
{
    const cyclicACMIPolyPatch& cpp =
        static_cast<const cyclicACMIPolyPatch&>(mesh_.boundaryMesh()[patch()]);

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
        couple = cpp.pointAMIAndFace(localFacei, direction, pos).first() >= 0;
        nonOverlap = !couple;
    }

    if (couple)
    {
        hitCyclicAMIPatch(cloud, td, direction);
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
void Foam::particle::hitCyclicRepeatAMIPatch
(
    TrackCloudType& cloud,
    trackingData& td,
    const vector& direction
)
{
    hitCyclicAMIPatch(cloud, td, direction);
}


template<class TrackCloudType>
void Foam::particle::hitProcessorPatch(TrackCloudType&, trackingData&)
{}


template<class TrackCloudType>
void Foam::particle::hitWallPatch(TrackCloudType&, trackingData&)
{}


// ************************************************************************* //
