/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::particle::stationaryTetGeometry
(
    vector& centre,
    vector& base,
    vector& vertex1,
    vector& vertex2
) const
{
    const triFace triIs(currentTetIndices().faceTriIs(mesh_));
    const vectorField& ccs = mesh_.cellCentres();
    const pointField& pts = mesh_.points();

    centre = ccs[celli_];
    base = pts[triIs[0]];
    vertex1 = pts[triIs[1]];
    vertex2 = pts[triIs[2]];
}


inline Foam::barycentricTensor Foam::particle::stationaryTetTransform() const
{
    vector centre, base, vertex1, vertex2;
    stationaryTetGeometry(centre, base, vertex1, vertex2);

    return barycentricTensor(centre, base, vertex1, vertex2);
}


inline void Foam::particle::movingTetGeometry
(
    const scalar fraction,
    Pair<vector>& centre,
    Pair<vector>& base,
    Pair<vector>& vertex1,
    Pair<vector>& vertex2
) const
{
    const triFace triIs(currentTetIndices().faceTriIs(mesh_));
    const pointField& ptsOld = mesh_.oldPoints();
    const pointField& ptsNew = mesh_.points();

    // !!! <-- We would be better off using mesh_.cellCentres() here. However,
    // we need to put a mesh_.oldCellCentres() method in for this to work. The
    // values obtained from the mesh and those obtained from the cell do not
    // necessarily match. See mantis #1993.
    const vector ccOld = mesh_.cells()[celli_].centre(ptsOld, mesh_.faces());
    const vector ccNew = mesh_.cells()[celli_].centre(ptsNew, mesh_.faces());

    // Old and new points and cell centres are not sub-cycled. If we are sub-
    // cycling, then we have to account for the timestep change here by
    // modifying the fractions that we take of the old and new geometry.
    const Pair<scalar> s = stepFractionSpan();
    const scalar f0 = s[0] + stepFraction_*s[1], f1 = fraction*s[1];

    centre[0] = ccOld + f0*(ccNew - ccOld);
    base[0] = ptsOld[triIs[0]] + f0*(ptsNew[triIs[0]] - ptsOld[triIs[0]]);
    vertex1[0] = ptsOld[triIs[1]] + f0*(ptsNew[triIs[1]] - ptsOld[triIs[1]]);
    vertex2[0] = ptsOld[triIs[2]] + f0*(ptsNew[triIs[2]] - ptsOld[triIs[2]]);

    centre[1] = f1*(ccNew - ccOld);
    base[1] = f1*(ptsNew[triIs[0]] - ptsOld[triIs[0]]);
    vertex1[1] = f1*(ptsNew[triIs[1]] - ptsOld[triIs[1]]);
    vertex2[1] = f1*(ptsNew[triIs[2]] - ptsOld[triIs[2]]);
}


inline Foam::Pair<Foam::barycentricTensor> Foam::particle::movingTetTransform
(
    const scalar fraction
) const
{
    Pair<vector> centre, base, vertex1, vertex2;
    movingTetGeometry(fraction, centre, base, vertex1, vertex2);

    return
        Pair<barycentricTensor>
        (
            barycentricTensor(centre[0], base[0], vertex1[0], vertex2[0]),
            barycentricTensor(centre[1], base[1], vertex1[1], vertex2[1])
        );
}


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


void Foam::particle::patchData(vector& n, vector& U) const
{
    if (!onBoundaryFace())
    {
        FatalErrorInFunction
            << "Patch data was requested for a particle that isn't on a patch"
            << exit(FatalError);
    }

    if (mesh_.moving())
    {
        Pair<vector> centre, base, vertex1, vertex2;
        movingTetGeometry(1, centre, base, vertex1, vertex2);

        n = triPointRef(base[0], vertex1[0], vertex2[0]).normal();

        // Interpolate the motion of the three face vertices to the current
        // coordinates
        U =
            coordinates_.b()*base[1]
          + coordinates_.c()*vertex1[1]
          + coordinates_.d()*vertex2[1];

        // The movingTetGeometry method gives the motion as a displacement
        // across the time-step, so we divide by the time-step to get velocity
        U /= mesh_.time().deltaTValue();
    }
    else
    {
        vector centre, base, vertex1, vertex2;
        stationaryTetGeometry(centre, base, vertex1, vertex2);

        n = triPointRef(base, vertex1, vertex2).normal();

        U = Zero;
    }
}


// ************************************************************************* //
