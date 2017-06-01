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

#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::particle::getNewParticleID() const
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


inline const Foam::polyMesh& Foam::particle::mesh() const
{
    return mesh_;
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


inline Foam::tetIndices Foam::particle::currentTetIndices() const
{
    return tetIndices(celli_, tetFacei_, tetPti_);
}


inline Foam::vector Foam::particle::normal() const
{
    return currentTetIndices().faceTri(mesh_).normal();
}


inline Foam::vector Foam::particle::oldNormal() const
{
    return currentTetIndices().oldFaceTri(mesh_).normal();
}


inline Foam::label Foam::particle::face() const
{
    return facei_;
}


inline bool Foam::particle::onFace() const
{
    return facei_ >= 0;
}


inline bool Foam::particle::onInternalFace() const
{
    return onFace() && mesh_.isInternalFace(facei_);
}


inline bool Foam::particle::onBoundaryFace() const
{
    return onFace() && !mesh_.isInternalFace(facei_);
}


inline Foam::vector Foam::particle::position() const
{
    if (mesh_.moving())
    {
        return movingTetTransform(0)[0] & coordinates_;
    }
    else
    {
        return tetTransform() & coordinates_;
    }
}


inline Foam::scalar& Foam::particle::stepFraction()
{
    return stepFraction_;
}


inline Foam::scalar Foam::particle::stepFraction() const
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


inline Foam::scalar Foam::particle::currentTime() const
{
    return mesh_.time().value() + stepFraction_*mesh_.time().deltaTValue();
}


inline Foam::label Foam::particle::patch(const label facei) const
{
    return mesh_.boundaryMesh().whichPatch(facei);
}


inline Foam::label Foam::particle::patchFace
(
    const label patchi,
    const label facei
) const
{
    return mesh_.boundaryMesh()[patchi].whichFace(facei);
}


// ************************************************************************* //
