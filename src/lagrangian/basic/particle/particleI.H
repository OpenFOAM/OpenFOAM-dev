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
#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::particle::getNewParticleIndex() const
{
    label id = particleCount_++;

    if (id == labelMax)
    {
        WarningInFunction
            << "Particle counter has overflowed. This might cause problems"
            << " when reconstructing particle tracks." << endl;
    }
    return id;
}


inline const Foam::barycentric& Foam::particle::coordinates() const
{
    return coordinates_;
}


inline Foam::label Foam::particle::cell() const
{
    return celli_;
}


inline Foam::label Foam::particle::tetFace() const
{
    return tetFacei_;
}


inline Foam::label Foam::particle::tetPt() const
{
    return tetPti_;
}


inline Foam::label Foam::particle::face() const
{
    return facei_;
}


inline Foam::label& Foam::particle::face()
{
    return facei_;
}


inline Foam::scalar Foam::particle::stepFraction() const
{
    return stepFraction_;
}


inline Foam::scalar& Foam::particle::stepFraction()
{
    return stepFraction_;
}


inline Foam::label Foam::particle::origProc() const
{
    return origProc_;
}


inline Foam::label& Foam::particle::origProc()
{
    return origProc_;
}


inline Foam::label Foam::particle::origId() const
{
    return origId_;
}


inline Foam::label& Foam::particle::origId()
{
    return origId_;
}


inline Foam::tetIndices Foam::particle::currentTetIndices
(
    const polyMesh& mesh
) const
{
    return tetIndices(celli_, tetFacei_, tetPti_);
}


inline Foam::vector Foam::particle::normal(const polyMesh& mesh) const
{
    return currentTetIndices(mesh).faceTri(mesh).normal();
}


inline bool Foam::particle::onFace() const
{
    return facei_ >= 0;
}


inline bool Foam::particle::onInternalFace(const polyMesh& mesh) const
{
    return onFace() && mesh.isInternalFace(facei_);
}


inline bool Foam::particle::onBoundaryFace(const polyMesh& mesh) const
{
    return onFace() && !mesh.isInternalFace(facei_);
}


inline Foam::label Foam::particle::patch(const polyMesh& mesh) const
{
    return onFace() ? mesh.boundaryMesh().whichPatch(facei_) : -1;
}


inline Foam::vector Foam::particle::position(const polyMesh& mesh) const
{
    return
        tracking::position
        (
            mesh,
            coordinates_,
            celli_,
            tetFacei_,
            tetPti_,
            stepFraction_
        );
}


inline void Foam::particle::reset(const scalar stepFraction)
{
    stepFraction_ = stepFraction;
    stepFractionBehind_ = 0;
    nTracksBehind_ = 0;
}


void Foam::particle::patchData
(
    const polyMesh& mesh,
    vector& normal,
    vector& displacement
) const
{
    if (!onBoundaryFace(mesh))
    {
        FatalErrorInFunction
            << "Patch data was requested for a particle that isn't on a patch"
            << exit(FatalError);
    }

    const Pair<vector> normalAndDisplacement =
        tracking::faceNormalAndDisplacement
        (
            mesh,
            coordinates_,
            celli_,
            tetFacei_,
            tetPti_,
            stepFraction_
        );

    normal = normalAndDisplacement.first();
    displacement = normalAndDisplacement.second();
}


// ************************************************************************* //
