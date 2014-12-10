/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "wallBoundedParticle.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class TrackData>
void Foam::wallBoundedParticle::patchInteraction
(
    TrackData& td,
    const scalar trackFraction
)
{
//    typedef TrackData::CloudType cloudType;
    typedef typename TrackData::cloudType::particleType particleType;

    particleType& p = static_cast<particleType&>(*this);
    p.hitFace(td);

    if (!internalFace(faceI_))
    {
        label origFaceI = faceI_;
        label patchI = patch(faceI_);

        // No action taken for tetPtI_ for tetFaceI_ here, handled by
        // patch interaction call or later during processor transfer.


        // Dummy tet indices. What to do here?
        tetIndices faceHitTetIs;

        if
        (
            !p.hitPatch
            (
                mesh_.boundaryMesh()[patchI],
                td,
                patchI,
                trackFraction,
                faceHitTetIs
            )
        )
        {
            // Did patch interaction model switch patches?
            // Note: recalculate meshEdgeStart_, diagEdge_!
            if (faceI_ != origFaceI)
            {
                patchI = patch(faceI_);
            }

            const polyPatch& patch = mesh_.boundaryMesh()[patchI];

            if (isA<wedgePolyPatch>(patch))
            {
                p.hitWedgePatch
                (
                    static_cast<const wedgePolyPatch&>(patch), td
                );
            }
            else if (isA<symmetryPlanePolyPatch>(patch))
            {
                p.hitSymmetryPlanePatch
                (
                    static_cast<const symmetryPlanePolyPatch&>(patch), td
                );
            }
            else if (isA<symmetryPolyPatch>(patch))
            {
                p.hitSymmetryPatch
                (
                    static_cast<const symmetryPolyPatch&>(patch), td
                );
            }
            else if (isA<cyclicPolyPatch>(patch))
            {
                p.hitCyclicPatch
                (
                    static_cast<const cyclicPolyPatch&>(patch), td
                );
            }
            else if (isA<processorPolyPatch>(patch))
            {
                p.hitProcessorPatch
                (
                    static_cast<const processorPolyPatch&>(patch), td
                );
            }
            else if (isA<wallPolyPatch>(patch))
            {
                p.hitWallPatch
                (
                    static_cast<const wallPolyPatch&>(patch), td, faceHitTetIs
                );
            }
            else
            {
                p.hitPatch(patch, td);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Track particle to a given position and returns 1.0 if the
//  trajectory is completed without hitting a face otherwise
//  stops at the face and returns the fraction of the trajectory
//  completed.
//  on entry 'stepFraction()' should be set to the fraction of the
//  time-step at which the tracking starts.
template<class TrackData>
Foam::scalar Foam::wallBoundedParticle::trackToEdge
(
    TrackData& td,
    const vector& endPosition
)
{
    // Are we on a track face? If not we do a topological walk.

    // Particle:
    // - cell_              always set
    // - tetFace_, tetPt_   always set (these identify tet particle is in)
    // - optionally meshEdgeStart_ or  diagEdge_ set (edge particle is on)

    //checkInside();
    //checkOnTriangle(position());
    //if (meshEdgeStart_ != -1 || diagEdge_ != -1)
    //{
    //    checkOnEdge();
    //}

    scalar trackFraction = 0.0;

    if (!td.isWallPatch_[tetFace()])
    {
        // Don't track across face. Just walk in cell. Particle is on
        // mesh edge (as indicated by meshEdgeStart_).
        const edge meshEdge(currentEdge());

        // If internal face check whether to go to neighbour cell or just
        // check to the other internal tet on the edge.
        if (mesh_.isInternalFace(tetFace()))
        {
            label nbrCellI =
            (
                cellI_ == mesh_.faceOwner()[faceI_]
              ? mesh_.faceNeighbour()[faceI_]
              : mesh_.faceOwner()[faceI_]
            );
            // Check angle to nbrCell tet. Is it in the direction of the
            // endposition? I.e. since volume of nbr tet is positive the
            // tracking direction should be into the tet.
            tetIndices nbrTi(nbrCellI, tetFaceI_, tetPtI_, mesh_);
            if ((nbrTi.faceTri(mesh_).normal() & (endPosition-position())) < 0)
            {
                // Change into nbrCell. No need to change tetFace, tetPt.
                //Pout<< "    crossed from cell:" << cellI_
                //    << " into " << nbrCellI << endl;
                cellI_ = nbrCellI;
                patchInteraction(td, trackFraction);
            }
            else
            {
                // Walk to other face on edge. Changes tetFace, tetPt but not
                // cell.
                crossEdgeConnectedFace(meshEdge);
                patchInteraction(td, trackFraction);
            }
        }
        else
        {
            // Walk to other face on edge. This might give loop since
            // particle should have been removed?
            crossEdgeConnectedFace(meshEdge);
            patchInteraction(td, trackFraction);
        }
    }
    else
    {
        // We're inside a tet on the wall. Check if the current tet is
        // the one to cross. If not we cross into the neighbouring triangle.

        if (mesh_.isInternalFace(tetFace()))
        {
            FatalErrorIn
            (
                "wallBoundedParticle::trackToEdge"
                "(TrackData&, const vector&)"
            )   << "Can only track on boundary faces."
                << " Face:" << tetFace()
                << " at:" << mesh_.faceCentres()[tetFace()]
                << abort(FatalError);
        }

        point projectedEndPosition = endPosition;
        // Remove normal component
        {
            const triFace tri(currentTetIndices().faceTriIs(mesh_));
            vector n = tri.normal(mesh_.points());
            n /= mag(n);
            const point& basePt = mesh_.points()[tri[0]];
            projectedEndPosition -= ((projectedEndPosition-basePt)&n)*n;
        }


        bool doTrack = false;
        if (meshEdgeStart_ == -1 && diagEdge_ == -1)
        {
            // We're starting and not yet on an edge.
            doTrack = true;
        }
        else
        {
            // See if the current triangle has got a point on the
            // correct side of the edge.
            doTrack = isTriAlongTrack(projectedEndPosition);
        }


        if (doTrack)
        {
            // Track across triangle. Return triangle edge crossed.
            label triEdgeI = -1;
            trackFraction = trackFaceTri(projectedEndPosition, triEdgeI);

            if (triEdgeI == -1)
            {
                // Reached endpoint
                //checkInside();
                diagEdge_ = -1;
                meshEdgeStart_ = -1;
                return trackFraction;
            }

            const tetIndices ti(currentTetIndices());

            // Triangle (faceTriIs) gets constructed from
            //    f[faceBasePtI_],
            //    f[facePtAI_],
            //    f[facePtBI_]
            //
            // So edge indices are:
            // 0 : edge between faceBasePtI_ and facePtAI_
            // 1 : edge between facePtAI_ and facePtBI_ (is always a real edge)
            // 2 : edge between facePtBI_ and faceBasePtI_

            const Foam::face& f = mesh_.faces()[ti.face()];
            const label fp0 = ti.faceBasePt();

            if (triEdgeI == 0)
            {
                if (ti.facePtA() == f.fcIndex(fp0))
                {
                    //Pout<< "Real edge." << endl;
                    diagEdge_ = -1;
                    meshEdgeStart_ = fp0;
                    //checkOnEdge();
                    crossEdgeConnectedFace(currentEdge());
                    patchInteraction(td, trackFraction);
                }
                else if (ti.facePtA() == f.rcIndex(fp0))
                {
                    //Note: should not happen since boundary face so owner
                    //Pout<< "Real edge." << endl;
                    FatalErrorIn("shold not happend") << info()
                        << abort(FatalError);

                    diagEdge_ = -1;
                    meshEdgeStart_ = f.rcIndex(fp0);
                    //checkOnEdge();
                    crossEdgeConnectedFace(currentEdge());
                    patchInteraction(td, trackFraction);
                }
                else
                {
                    // Get index of triangle on other side of edge.
                    diagEdge_ = ti.facePtA()-fp0;
                    if (diagEdge_ < 0)
                    {
                        diagEdge_ += f.size();
                    }
                    meshEdgeStart_ = -1;
                    //checkOnEdge();
                    crossDiagonalEdge();
                }
            }
            else if (triEdgeI == 1)
            {
                //Pout<< "Real edge." << endl;
                diagEdge_ = -1;
                meshEdgeStart_ = ti.facePtA();
                //checkOnEdge();
                crossEdgeConnectedFace(currentEdge());
                patchInteraction(td, trackFraction);
            }
            else // if (triEdgeI == 2)
            {
                if (ti.facePtB() == f.rcIndex(fp0))
                {
                    //Pout<< "Real edge." << endl;
                    diagEdge_ = -1;
                    meshEdgeStart_ = ti.facePtB();
                    //checkOnEdge();
                    crossEdgeConnectedFace(currentEdge());
                    patchInteraction(td, trackFraction);
                }
                else if (ti.facePtB() == f.fcIndex(fp0))
                {
                    //Note: should not happen since boundary face so owner
                    //Pout<< "Real edge." << endl;
                    FatalErrorIn("shold not happend") << info()
                        << abort(FatalError);

                    diagEdge_ = -1;
                    meshEdgeStart_ = fp0;
                    //checkOnEdge();
                    crossEdgeConnectedFace(currentEdge());
                    patchInteraction(td, trackFraction);
                }
                else
                {
                    //Pout<< "Triangle edge." << endl;
                    // Get index of triangle on other side of edge.
                    diagEdge_ = ti.facePtB()-fp0;
                    if (diagEdge_ < 0)
                    {
                        diagEdge_ += f.size();
                    }
                    meshEdgeStart_ = -1;
                    //checkOnEdge();
                    crossDiagonalEdge();
                }
            }
        }
        else
        {
            // Current tet is not the right one. Check the neighbour tet.

            if (meshEdgeStart_ != -1)
            {
                // Particle is on mesh edge so change into other face on cell
                crossEdgeConnectedFace(currentEdge());
                //checkOnEdge();
                patchInteraction(td, trackFraction);
            }
            else
            {
                // Particle is on diagonal edge so change into the other
                // triangle.
                crossDiagonalEdge();
                //checkOnEdge();
            }
        }
    }

    //checkInside();

    return trackFraction;
}


template<class TrackData>
bool Foam::wallBoundedParticle::hitPatch
(
    const polyPatch&,
    TrackData& td,
    const label patchI,
    const scalar trackFraction,
    const tetIndices& tetIs
)
{
    // Disable generic patch interaction
    return false;
}


template<class TrackData>
void Foam::wallBoundedParticle::hitWedgePatch
(
    const wedgePolyPatch& pp,
    TrackData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


template<class TrackData>
void Foam::wallBoundedParticle::hitSymmetryPlanePatch
(
    const symmetryPlanePolyPatch& pp,
    TrackData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


template<class TrackData>
void Foam::wallBoundedParticle::hitSymmetryPatch
(
    const symmetryPolyPatch& pp,
    TrackData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


template<class TrackData>
void Foam::wallBoundedParticle::hitCyclicPatch
(
    const cyclicPolyPatch& pp,
    TrackData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


template<class TrackData>
void Foam::wallBoundedParticle::hitProcessorPatch
(
    const processorPolyPatch& pp,
    TrackData& td
)
{
    // Switch particle
    td.switchProcessor = true;

    // Adapt edgeStart_ for other side.
    // E.g. if edgeStart_ is 1 then the edge is between vertex 1 and 2 so
    // on the other side between 2 and 3 so edgeStart_ should be
    // f.size()-edgeStart_-1.

    const Foam::face& f = mesh_.faces()[face()];

    if (meshEdgeStart_ != -1)
    {
        meshEdgeStart_ = f.size()-meshEdgeStart_-1;
    }
    else
    {
        // diagEdge_ is relative to faceBasePt
        diagEdge_ = f.size()-diagEdge_;
    }
}


template<class TrackData>
void Foam::wallBoundedParticle::hitWallPatch
(
    const wallPolyPatch& wpp,
    TrackData& td,
    const tetIndices&
)
{}


template<class TrackData>
void Foam::wallBoundedParticle::hitPatch
(
    const polyPatch& wpp,
    TrackData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


template<class CloudType>
void Foam::wallBoundedParticle::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);

    IOField<label> meshEdgeStart
    (
        c.fieldIOobject("meshEdgeStart", IOobject::MUST_READ)
    );

    IOField<label> diagEdge
    (
        c.fieldIOobject("diagEdge_", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, diagEdge);

    label i = 0;
    forAllIter(typename CloudType, c, iter)
    {
        iter().meshEdgeStart_ = meshEdgeStart[i];
        iter().diagEdge_ = diagEdge[i];
        i++;
    }
}


template<class CloudType>
void Foam::wallBoundedParticle::writeFields(const CloudType& c)
{
    particle::writeFields(c);

    label np =  c.size();

    IOField<label> meshEdgeStart
    (
        c.fieldIOobject("meshEdgeStart", IOobject::NO_READ),
        np
    );
    IOField<label> diagEdge
    (
        c.fieldIOobject("diagEdge", IOobject::NO_READ),
        np
    );

    label i = 0;
    forAllConstIter(typename CloudType, c, iter)
    {
        meshEdgeStart[i] = iter().meshEdgeStart_;
        diagEdge[i] = iter().diagEdge_;
        i++;
    }

    meshEdgeStart.write();
    diagEdge.write();
}


// ************************************************************************* //
