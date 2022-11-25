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
#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::particle::stationaryTetGeometry
(
    const polyMesh& mesh,
    vector& centre,
    vector& base,
    vector& vertex1,
    vector& vertex2
) const
{
    const triFace triIs(currentTetIndices(mesh).faceTriIs(mesh));
    const vectorField& ccs = mesh.cellCentres();
    const pointField& pts = mesh.points();

    centre = ccs[celli_];
    base = pts[triIs[0]];
    vertex1 = pts[triIs[1]];
    vertex2 = pts[triIs[2]];
}


inline Foam::barycentricTensor Foam::particle::stationaryTetTransform
(
    const polyMesh& mesh
) const
{
    vector centre, base, vertex1, vertex2;
    stationaryTetGeometry(mesh, centre, base, vertex1, vertex2);

    return barycentricTensor(centre, base, vertex1, vertex2);
}


inline void Foam::particle::movingTetGeometry
(
    const polyMesh& mesh,
    const scalar fraction,
    Pair<vector>& centre,
    Pair<vector>& base,
    Pair<vector>& vertex1,
    Pair<vector>& vertex2
) const
{
    const triFace triIs(currentTetIndices(mesh).faceTriIs(mesh));
    const pointField& ptsOld = mesh.oldPoints();
    const pointField& ptsNew = mesh.points();
    const vector ccOld = mesh.oldCellCentres()[celli_];
    const vector ccNew = mesh.cellCentres()[celli_];

    const scalar f0 = stepFraction_, f1 = fraction;

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
    const polyMesh& mesh,
    const scalar fraction
) const
{
    Pair<vector> centre, base, vertex1, vertex2;
    movingTetGeometry(mesh, fraction, centre, base, vertex1, vertex2);

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


inline Foam::barycentricTensor Foam::particle::currentTetTransform
(
    const polyMesh& mesh
) const
{
    if (mesh.moving() && stepFraction_ != 1)
    {
        return movingTetTransform(mesh, 0)[0];
    }
    else
    {
        return stationaryTetTransform(mesh);
    }
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
    return currentTetTransform(mesh) & coordinates_;
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

    if (mesh.moving() && stepFraction_ != 1)
    {
        Pair<vector> centre, base, vertex1, vertex2;
        movingTetGeometry(mesh, 1, centre, base, vertex1, vertex2);

        normal = triPointRef(base[0], vertex1[0], vertex2[0]).normal();

        // Interpolate the motion of the three face vertices to the current
        // coordinates
        displacement =
            coordinates_.b()*base[1]
          + coordinates_.c()*vertex1[1]
          + coordinates_.d()*vertex2[1];
    }
    else
    {
        vector centre, base, vertex1, vertex2;
        stationaryTetGeometry(mesh, centre, base, vertex1, vertex2);

        normal = triPointRef(base, vertex1, vertex2).normal();

        displacement = Zero;
    }
}


// ************************************************************************* //
