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
void Foam::particle::prepareForParallelTransfer
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    if (td.sendFromPatch == patch(td.mesh))
    {
        prepareForProcessorTransfer(td);
    }
    else
    {
        prepareForNonConformalCyclicTransfer
        (
            td.mesh,
            td.sendFromPatch,
            td.sendToPatchFace
        );
    }
}


template<class TrackCloudType>
void Foam::particle::correctAfterParallelTransfer
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    const polyPatch& pp = td.mesh.boundaryMesh()[td.sendToPatch];

    if (isA<processorPolyPatch>(pp))
    {
        correctAfterProcessorTransfer(td);
    }
    else if (isA<nonConformalCyclicPolyPatch>(pp))
    {
        correctAfterNonConformalCyclicTransfer(td.mesh, td.sendToPatch);
    }
    else
    {
        FatalErrorInFunction
            << "Transfer patch type not recognised"
            << exit(FatalError);
    }
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
    else if (onInternalFace(td.mesh))
    {
        changeCell(td.mesh);
    }
    else if (onBoundaryFace(td.mesh))
    {
        forAll(cloud.patchNonConformalCyclicPatches()[p.patch(td.mesh)], i)
        {
            if
            (
                p.hitNonConformalCyclicPatch
                (
                    displacement,
                    fraction,
                    cloud.patchNonConformalCyclicPatches()[p.patch(td.mesh)][i],
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
            const polyPatch& patch = td.mesh.boundaryMesh()[p.patch(td.mesh)];

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
                p.hitBasicPatch(cloud, ttd);
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

    const scalar f = trackToFace(td.mesh, displacement, fraction);

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
void Foam::particle::hitSymmetryPatch(TrackCloudType&, trackingData& td)
{
    const vector nf = normal(td.mesh);
    transformProperties(transformer::rotation(I - 2.0*nf*nf));
}


template<class TrackCloudType>
void Foam::particle::hitCyclicPatch(TrackCloudType&, trackingData& td)
{
    const cyclicPolyPatch& cpp =
        static_cast<const cyclicPolyPatch&>
        (
            td.mesh.boundaryMesh()[patch(td.mesh)]
        );
    const cyclicPolyPatch& receiveCpp = cpp.nbrPatch();

    // Set the topology
    facei_ = tetFacei_ = cpp.transformGlobalFace(facei_);
    celli_ = td.mesh.faceOwner()[facei_];

    // See note in correctAfterParallelTransfer for tetPti addressing ...
    tetPti_ = td.mesh.faces()[tetFacei_].size() - 1 - tetPti_;

    // Reflect to account for the change of triangle orientation in the new cell
    reflect();

    // Transform the properties
    if (receiveCpp.transform().transformsPosition())
    {
        transformProperties(receiveCpp.transform());
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
        (td.mesh.boundaryMesh()[patchi]);

    const point sendPos = position(td.mesh);

    // Get the send patch data
    vector sendNormal, sendDisplacement;
    patchData(td.mesh, sendNormal, sendDisplacement);

    // Project the particle through the non-conformal patch
    point receivePos;
    const remote receiveProcFace =
        nccpp.ray
        (
            stepFraction_,
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
        td.sendToPatchFace = receiveProcFace.elementi;

        return true;
    }

    // If both sides are on the same process, then do the local transfer
    prepareForNonConformalCyclicTransfer
    (
        td.mesh,
        nccpp.index(),
        receiveProcFace.elementi
    );
    correctAfterNonConformalCyclicTransfer
    (
        td.mesh,
        nccpp.nbrPatchID()
    );

    return true;
}


template<class TrackCloudType>
void Foam::particle::hitProcessorPatch(TrackCloudType& cloud, trackingData& td)
{
    td.sendToProc = cloud.patchNbrProc()[patch(td.mesh)];
    td.sendFromPatch = patch(td.mesh);
    td.sendToPatch = cloud.patchNbrProcPatch()[patch(td.mesh)];
    td.sendToPatchFace =
        td.mesh.boundaryMesh()[patch(td.mesh)].whichFace(face());
}


template<class TrackCloudType>
void Foam::particle::hitWallPatch(TrackCloudType&, trackingData&)
{}


template<class TrackCloudType>
void Foam::particle::hitBasicPatch(TrackCloudType&, trackingData& td)
{
    td.keepParticle = false;
}


// ************************************************************************* //
