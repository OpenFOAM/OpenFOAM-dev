/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
#include "tracking.H"
#include "polyTopoChangeMap.H"
#include "transform.H"
#include "treeDataCell.H"
#include "indexedOctree.H"
#include "cubicEqn.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::particle::particleCount_ = 0;

namespace Foam
{
    defineTypeNameAndDebug(particle, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particle::particle
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label facei
)
:
    coordinates_(coordinates),
    celli_(celli),
    tetFacei_(tetFacei),
    tetPti_(tetPti),
    facei_(facei),
    stepFraction_(1),
    stepFractionBehind_(0),
    nTracksBehind_(0),
    origProc_(Pstream::myProcNo()),
    origId_(getNewParticleIndex())
{}


Foam::particle::particle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    label& nLocateBoundaryHits
)
:
    coordinates_(- vGreat, - vGreat, - vGreat, - vGreat),
    celli_(celli),
    tetFacei_(-1),
    tetPti_(-1),
    facei_(-1),
    stepFraction_(1),
    stepFractionBehind_(0),
    nTracksBehind_(0),
    origProc_(Pstream::myProcNo()),
    origId_(getNewParticleIndex())
{
    if (!locate(mesh, position, celli))
    {
        nLocateBoundaryHits ++;
    }
}


Foam::particle::particle(const particle& p)
:
    coordinates_(p.coordinates_),
    celli_(p.celli_),
    tetFacei_(p.tetFacei_),
    tetPti_(p.tetPti_),
    facei_(p.facei_),
    stepFraction_(p.stepFraction_),
    stepFractionBehind_(p.stepFractionBehind_),
    nTracksBehind_(p.nTracksBehind_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::particle::locate
(
    const polyMesh& mesh,
    const vector& position,
    label celli
)
{
    celli_ = celli;

    return
        tracking::locate
        (
            mesh, position,
            coordinates_, celli_, tetFacei_, tetPti_, 1,
            debug
          ? static_cast<const string&>(string("Particle " + name(origId())))
          : NullObjectRef<string>()
        );
}


Foam::scalar Foam::particle::track
(
    const polyMesh& mesh,
    const vector& displacement,
    const scalar fraction
)
{
    const Tuple2<bool, scalar> onBoundaryAndF =
        tracking::toBoundary
        (
            mesh, displacement, fraction,
            coordinates_, celli_, tetFacei_, tetPti_, stepFraction_,
            stepFractionBehind_, nTracksBehind_,
            debug
          ? static_cast<const string&>(string("Particle " + name(origId())))
          : NullObjectRef<string>()
        );

    facei_ = onBoundaryAndF.first() ? tetFacei_ : -1;

    return onBoundaryAndF.second();
}


Foam::scalar Foam::particle::trackToCell
(
    const polyMesh& mesh,
    const vector& displacement,
    const scalar fraction
)
{
    const Tuple2<bool, scalar> inNextCellAndF =
        tracking::toCell
        (
            mesh, displacement, fraction,
            coordinates_, celli_, tetFacei_, tetPti_, stepFraction_,
            stepFractionBehind_, nTracksBehind_,
            debug
          ? static_cast<const string&>(string("Particle " + name(origId())))
          : NullObjectRef<string>()
        );

    facei_ = inNextCellAndF.first() ? tetFacei_ : -1;

    return inNextCellAndF.second();
}


Foam::scalar Foam::particle::trackToFace
(
    const polyMesh& mesh,
    const vector& displacement,
    const scalar fraction
)
{
    const Tuple2<bool, scalar> onFaceAndF =
        tracking::toFace
        (
            mesh, displacement, fraction,
            coordinates_, celli_, tetFacei_, tetPti_, stepFraction_,
            stepFractionBehind_, nTracksBehind_,
            debug
          ? static_cast<const string&>(string("Particle " + name(origId())))
          : NullObjectRef<string>()
        );

    facei_ = onFaceAndF.first() ? tetFacei_ : -1;

    return onFaceAndF.second();
}


Foam::vector Foam::particle::deviationFromMeshCentre
(
    const polyMesh& mesh
) const
{
    if (cmptMin(mesh.geometricD()) == -1)
    {
        vector pos = position(mesh), posC = pos;
        meshTools::constrainToMeshCentre(mesh, posC);
        return pos - posC;
    }
    else
    {
        return vector::zero;
    }
}


void Foam::particle::transformProperties(const transformer&)
{}


void Foam::particle::prepareForProcessorTransfer(trackingData& td)
{
    const processorPolyPatch& ppp =
        refCast<const processorPolyPatch>
        (
            td.mesh.boundaryMesh()[td.sendFromPatch]
        );

    tracking::inProcessor(ppp, celli_, tetFacei_);
}


void Foam::particle::correctAfterProcessorTransfer(trackingData& td)
{
    const processorPolyPatch& ppp =
        refCast<const processorPolyPatch>
        (
            td.mesh.boundaryMesh()[td.sendToPatch]
        );

    tracking::outProcessor(ppp, coordinates_, celli_, tetFacei_, tetPti_);

    // Remain on the face
    facei_ = tetFacei_;

    // Transform the properties
    if (ppp.transform().transformsPosition())
    {
        transformProperties(ppp.transform());
    }
}


void Foam::particle::prepareForNonConformalCyclicTransfer
(
    const polyMesh& mesh,
    const label sendFromPatch,
    const label sendToPatchFace,
    const vector& sendToPosition
)
{
    const nonConformalCyclicPolyPatch& nccpp =
        static_cast<const nonConformalCyclicPolyPatch&>
        (
            mesh.boundaryMesh()[sendFromPatch]
        );

    // Store the position in the barycentric data
    coordinates_ =
        barycentric
        (
            1 - cmptSum(sendToPosition),
            sendToPosition.x(),
            sendToPosition.y(),
            sendToPosition.z()
        );

    // Break the topology
    celli_ = -1;
    tetFacei_ = -1;
    tetPti_ = -1;

    // Store the local patch face in the face index
    facei_ = sendToPatchFace;

    // Transform the properties
    if (nccpp.transform().transformsPosition())
    {
        transformProperties(nccpp.nbrPatch().transform());
    }
}


void Foam::particle::correctAfterNonConformalCyclicTransfer
(
    const polyMesh& mesh,
    const label sendToPatch,
    labelList& patchNLocateBoundaryHits
)
{
    const nonConformalCyclicPolyPatch& nccpp =
        static_cast<const nonConformalCyclicPolyPatch&>
        (
            mesh.boundaryMesh()[sendToPatch]
        );

    // Get the position from the barycentric data
    const vector receivePos
    (
        coordinates_.b(),
        coordinates_.c(),
        coordinates_.d()
    );

    // Locate the particle on the receiving side
    const label celli = mesh.faceOwner()[facei_ + nccpp.origPatch().start()];
    if (!locate(mesh, receivePos, celli))
    {
        patchNLocateBoundaryHits[sendToPatch] ++;
    }

    // The particle must remain associated with a face for the tracking to
    // register as incomplete
    facei_ = tetFacei_;
}


void Foam::particle::prepareForInteractionListReferral
(
    const polyMesh& mesh,
    const transformer& transform
)
{
    // Get the transformed position
    const vector pos = transform.invTransformPosition(position(mesh));

    // Break the topology
    celli_ = -1;
    tetFacei_ = -1;
    tetPti_ = -1;

    // Store the position in the barycentric data
    coordinates_ = barycentric(1 - cmptSum(pos), pos.x(), pos.y(), pos.z());

    // Transform the properties
    if (transform.transformsPosition())
    {
        transformProperties(inv(transform));
    }
}


void Foam::particle::correctAfterInteractionListReferral
(
    const polyMesh& mesh,
    const label celli
)
{
    // Get the position from the barycentric data
    const vector pos(coordinates_.b(), coordinates_.c(), coordinates_.d());

    // Create some arbitrary topology for the supplied cell
    celli_ = celli;
    tetFacei_ = mesh.cells()[celli_][0];
    tetPti_ = 1;

    // Calculate the coordinates for the assumed topology. This isn't likely to
    // be correct; the particle is probably not in this tet. It will, however,
    // generate the correct vector when the position method is called. A
    // referred particle should never be tracked, so this approximate topology
    // is good enough. By using the nearby cell we minimise the error
    // associated with the incorrect topology.
    coordinates_ =
        tracking::coordinates
        (
            mesh, pos,
            celli_, tetFacei_, tetPti_, stepFraction_
        );
}


Foam::label Foam::particle::procTetPt
(
    const polyMesh& mesh,
    const polyMesh& procMesh,
    const label procCell,
    const label procTetFace
) const
{
    // The tet point on the procMesh differs from the current tet point if the
    // mesh and procMesh faces are of differing orientation. The change is the
    // same as in particle::correctAfterParallelTransfer.

    if
    (
        (mesh.faceOwner()[tetFacei_] == celli_)
     == (procMesh.faceOwner()[procTetFace] == procCell)
    )
    {
        return tetPti_;
    }
    else
    {
        return procMesh.faces()[procTetFace].size() - 1 - tetPti_;
    }
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //
//

bool Foam::operator==(const particle& pA, const particle& pB)
{
    return (pA.origProc() == pB.origProc() && pA.origId() == pB.origId());
}


bool Foam::operator!=(const particle& pA, const particle& pB)
{
    return !(pA == pB);
}


// ************************************************************************* //
