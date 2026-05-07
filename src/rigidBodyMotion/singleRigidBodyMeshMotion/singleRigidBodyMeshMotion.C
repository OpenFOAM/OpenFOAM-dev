/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2026 OpenFOAM Foundation
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

#include "singleRigidBodyMeshMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "pointDist.H"
#include "pointConstraints.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singleRigidBodyMeshMotion, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        singleRigidBodyMeshMotion,
        dictionary
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::singleRigidBodyMeshMotion::calcScale()
{
    const pointMesh& pMesh = pointMesh::New(mesh());

    const pointDist pDist(pMesh, patchSet_, points0(), do_);

    // One before the inner distance, zero after the outer distance, and a
    // linear variation in between
    scale_.primitiveFieldRef() =
        min
        (
            max
            (
                (do_ - pDist.primitiveField())/(do_ - di_),
                scalar(0)
            ),
            scalar(1)
        );

    // Convert the linear variation into a smooth cosine
    scale_.primitiveFieldRef() =
        min
        (
            max
            (
                0.5
              - 0.5
               *cos(scale_.primitiveField()
               *Foam::constant::mathematical::pi),
                scalar(0)
            ),
            scalar(1)
        );

    pointConstraints::New(pMesh).constrain(scale_);

    scale_.write();
}


Foam::tmp<Foam::pointField> Foam::singleRigidBodyMeshMotion::calcPoints
(
    const septernion& s,
    const vector& CofG
) const
{
    const pointField& points0 = this->points0();

    tmp<pointField> tpoints(new pointField(points0));
    pointField& points = tpoints.ref();

    forAll(points, pointi)
    {
        // Move non-stationary points
        if (scale_[pointi] > small)
        {
            // Use solid-body motion where scale = 1
            if (scale_[pointi] > 1 - small)
            {
                points[pointi] =
                    CofG + s.transformPoint(points0[pointi] - CofG);
            }
            // Slerp septernion interpolation
            else
            {
                const septernion ss(slerp(septernion::I, s, scale_[pointi]));

                points[pointi] =
                    CofG + ss.transformPoint(points0[pointi] - CofG);
            }
        }
    }

    return tpoints;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singleRigidBodyMeshMotion::singleRigidBodyMeshMotion
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    points0MotionSolver(name, mesh, dict, typeName),
    SBMFPtr_(solidBodyMotionFunction::New(dict, mesh.time())),
    patches_(wordReList(dict.lookup("patches"))),
    patchSet_(mesh.boundary().patchSet(patches_)),
    CofG_(dict.lookup("CofG")),
    di_(dict.lookup<scalar>("innerDistance")),
    do_(dict.lookup<scalar>("outerDistance")),
    scale_
    (
        IOobject
        (
            "motionScale",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, 0)
    )
{
    // Calculate scaling factor everywhere
    calcScale();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singleRigidBodyMeshMotion::~singleRigidBodyMeshMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::singleRigidBodyMeshMotion::newPoints()
{
    return calcPoints(SBMFPtr_().transformation(), CofG_);
}


void Foam::singleRigidBodyMeshMotion::topoChange
(
    const polyTopoChangeMap& map
)
{
    // Pending implementation of the inverse transformation of points0
    NotImplemented;
}


void Foam::singleRigidBodyMeshMotion::mapMesh(const polyMeshMap& map)
{
    // Reset the points0 to the current point locations
    points0MotionSolver::mapMesh(map);

    // Recalculate the scale. It will already have been resized by the mapping.
    calcScale();

    // Inverse-transform the points0-s so that the transformation still applies
    // correctly. This is fine because all the points have changed during the
    // mesh-to-mesh mapping. We are not concerned with returning exactly to the
    // original point locations, as those locations now relate to an entirely
    // different mesh. If we really want to return to those locations, then we
    // need to map back to the original mesh.
    points0_.primitiveFieldRef() =
        calcPoints(inv(SBMFPtr_().transformation()), -CofG_);

    // Marks points0 as to be (re-) written
    points0_.writeOpt() = IOobject::AUTO_WRITE;
    points0_.instance() = mesh().time().name();
}


// ************************************************************************* //
