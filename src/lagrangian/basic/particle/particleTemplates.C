/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "particle.H"
#include "IOPosition.H"

#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "nonConformalCyclicPolyPatch.H"
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

    typeIOobject<IOField<label>> procIO
    (
        c.fieldIOobject("origProcId", IOobject::MUST_READ)
    );

    bool haveFile = procIO.headerOk();

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
    const vector& displacement,
    const scalar fraction,
    TrackCloudType& cloud,
    trackingData& td
)
{
    if (debug)
    {
        Info << "Particle " << origId() << nl << FUNCTION_NAME << nl << endl;
    }

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
        forAll(cloud.patchNonConformalCyclicPatches()[p.patch()], i)
        {
            if
            (
                p.hitNonConformalCyclicPatch
                (
                    displacement,
                    fraction,
                    cloud.patchNonConformalCyclicPatches()[p.patch()][i],
                    cloud,
                    ttd
                )
            )
            {
                return;
            }
        }

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
            else if (isA<cyclicAMIPolyPatch>(patch))
            {
                p.hitCyclicAMIPatch(displacement, fraction, cloud, ttd);
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
Foam::scalar Foam::particle::trackToAndHitFace
(
    const vector& displacement,
    const scalar fraction,
    TrackCloudType& cloud,
    trackingData& td
)
{
    if (debug)
    {
        Info << "Particle " << origId() << nl << FUNCTION_NAME << nl << endl;
    }

    const scalar f = trackToFace(displacement, fraction);

    hitFace(displacement, fraction, cloud, td);

    return f;
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
    transformProperties(transformer::rotation(I - 2.0*nf*nf));
}


template<class TrackCloudType>
void Foam::particle::hitCyclicPatch(TrackCloudType&, trackingData&)
{
    const cyclicPolyPatch& cpp =
        static_cast<const cyclicPolyPatch&>(mesh_.boundaryMesh()[patch()]);
    const cyclicPolyPatch& receiveCpp = cpp.nbrPatch();

    // Set the topology
    facei_ = tetFacei_ = cpp.transformGlobalFace(facei_);
    celli_ = mesh_.faceOwner()[facei_];

    // See note in correctAfterParallelTransfer for tetPti addressing ...
    tetPti_ = mesh_.faces()[tetFacei_].size() - 1 - tetPti_;

    // Reflect to account for the change of triangle orientation in the new cell
    reflect();

    // Transform the properties
    if (receiveCpp.transform().transformsPosition())
    {
        transformProperties(receiveCpp.transform());
    }
}


template<class TrackCloudType>
void Foam::particle::hitCyclicAMIPatch
(
    const vector& displacement,
    const scalar fraction,
    TrackCloudType& cloud,
    trackingData& td
)
{
    const cyclicAMIPolyPatch& cpp =
        static_cast<const cyclicAMIPolyPatch&>(mesh_.boundaryMesh()[patch()]);
    const cyclicAMIPolyPatch& receiveCpp = cpp.nbrPatch();

    if (debug)
    {
        Info<< "Particle " << origId() << " crossing AMI from " << cpp.name()
            << " to " << receiveCpp.name() << endl << endl;
    }

    // Get the send patch data
    vector sendNormal, sendDisplacement;
    patchData(sendNormal, sendDisplacement);

    vector pos = position();

    const labelPair receiveIs =
        cpp.pointAMIAndFace
        (
            cpp.whichFace(facei_),
            displacement - fraction*sendDisplacement,
            pos
        );
    const label receiveAMIi = receiveIs.first();
    const label receiveFacei = receiveIs.second();

    // If the receiving face could not be found then issue a warning and remove
    // the particle
    if (receiveFacei < 0)
    {
        td.keepParticle = false;
        WarningInFunction
            << "Particle transfer from " << cyclicAMIPolyPatch::typeName
            << " patches " << cpp.name() << " to " << receiveCpp.name()
            << " failed at position " << pos << " and with displacement "
            << (displacement - fraction*sendDisplacement) << nl
            << "    A receiving face could not be found" << nl
            << "    The particle has been removed" << nl << endl;
        return;
    }

    // Set the topology
    facei_ = tetFacei_ = receiveFacei + receiveCpp.start();

    // Locate the particle on the receiving side
    locate
    (
        pos,
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
    vector displacementT = displacement;

    const transformer AMITransform =
        receiveCpp.owner()
      ? receiveCpp.AMITransforms()[receiveAMIi]
      : inv(cpp.AMITransforms()[receiveAMIi]);

    if (AMITransform.transformsPosition())
    {
        transformProperties(AMITransform);
        displacementT = AMITransform.transform(displacementT);
    }

    if (receiveCpp.transform().transformsPosition())
    {
        transformProperties(receiveCpp.transform());
        displacementT = receiveCpp.transform().transform(displacementT);
    }

    // If on a boundary and the displacement points into the receiving face
    // then issue a warning and remove the particle
    if (onBoundaryFace())
    {
        vector receiveNormal, receiveDisplacement;
        patchData(receiveNormal, receiveDisplacement);

        if (((displacementT - fraction*receiveDisplacement)&receiveNormal) > 0)
        {
            td.keepParticle = false;
            WarningInFunction
                << "Particle transfer from " << cyclicAMIPolyPatch::typeName
                << " patches " << cpp.name() << " to " << receiveCpp.name()
                << " failed at position " << pos << " and with displacement "
                << (displacementT - fraction*receiveDisplacement) << nl
                << "    The displacement points into both the source and "
                << "receiving faces, so the tracking cannot proceed" << nl
                << "    The particle has been removed" << nl << endl;
            return;
        }
    }
}


template<class TrackCloudType>
bool Foam::particle::hitNonConformalCyclicPatch
(
    const vector& displacement,
    const scalar fraction,
    const label patchi,
    TrackCloudType& cloud,
    trackingData& td
)
{
    const nonConformalCyclicPolyPatch& nccpp =
        static_cast<const nonConformalCyclicPolyPatch&>
        (mesh_.boundaryMesh()[patchi]);

    const point sendPos = this->position();

    // Get the send patch data
    vector sendNormal, sendDisplacement;
    patchData(sendNormal, sendDisplacement);

    // Project the particle through the non-conformal patch
    point receivePos;
    const patchToPatch::procFace receiveProcFace =
        nccpp.ray
        (
            stepFractionSpan()[0] + stepFraction_*stepFractionSpan()[1],
            nccpp.origPatch().whichFace(facei_),
            sendPos,
            displacement - fraction*sendDisplacement,
            receivePos
        );

    // If we didn't hit anything then this particle is assumed to project to
    // the orig patch, or another different non-conformal patch. Return, so
    // these can be tried.
    if (receiveProcFace.proci == -1) return false;

    // If we are transferring between processes then mark as such and return.
    // The cloud will handle all processor transfers as a single batch.
    if (receiveProcFace.proci != Pstream::myProcNo())
    {
        td.sendFromPatch = nccpp.index();
        td.sendToProc = receiveProcFace.proci;
        td.sendToPatch = nccpp.nbrPatchID();
        td.sendToPatchFace = receiveProcFace.facei;

        return true;
    }

    // If both sides are on the same process, then do the local transfer
    prepareForNonConformalCyclicTransfer
    (
        nccpp.index(),
        receiveProcFace.facei
    );
    correctAfterNonConformalCyclicTransfer
    (
        nccpp.nbrPatchID()
    );

    return true;
}


template<class TrackCloudType>
void Foam::particle::hitProcessorPatch(TrackCloudType& cloud, trackingData& td)
{
    td.sendToProc = cloud.patchNbrProc()[patch()];
    td.sendFromPatch = patch();
    td.sendToPatch = cloud.patchNbrProcPatch()[patch()];
    td.sendToPatchFace = mesh().boundaryMesh()[patch()].whichFace(face());
}


template<class TrackCloudType>
void Foam::particle::hitWallPatch(TrackCloudType&, trackingData&)
{}


// ************************************************************************* //
