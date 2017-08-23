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


inline Foam::label Foam::particle::face() const
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


inline Foam::Pair<Foam::scalar> Foam::particle::stepFractionSpan() const
{
    if (mesh_.time().subCycling())
    {
        const TimeState& tsNew = mesh_.time();
        const TimeState& tsOld = mesh_.time().prevTimeState();

        const scalar tFrac =
        (
            (tsNew.value() - tsNew.deltaTValue())
          - (tsOld.value() - tsOld.deltaTValue())
        )/tsOld.deltaTValue();

        const scalar dtFrac = tsNew.deltaTValue()/tsOld.deltaTValue();

        return Pair<scalar>(tFrac, dtFrac);
    }
    else
    {
        return Pair<scalar>(0, 1);
    }
}


inline Foam::scalar Foam::particle::currentTimeFraction() const
{
    const Pair<scalar> s = stepFractionSpan();

    return s[0] + stepFraction_*s[1];
}


inline Foam::tetIndices Foam::particle::currentTetIndices() const
{
    return tetIndices(celli_, tetFacei_, tetPti_);
}


inline Foam::barycentricTensor Foam::particle::currentTetTransform() const
{
    if (mesh_.moving())
    {
        return movingTetTransform(0)[0];
    }
    else
    {
        return stationaryTetTransform();
    }
}


inline Foam::vector Foam::particle::normal() const
{
    return currentTetIndices().faceTri(mesh_).normal();
}


inline Foam::vector Foam::particle::oldNormal() const
{
    return currentTetIndices().oldFaceTri(mesh_).normal();
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


inline Foam::label Foam::particle::patch() const
{
    return onFace() ? mesh_.boundaryMesh().whichPatch(facei_) : -1;
}


inline Foam::vector Foam::particle::position() const
{
    return currentTetTransform() & coordinates_;
}


// ************************************************************************* //
