/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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

#include "interpolatingSolidBodyMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interpolatingSolidBodyMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        interpolatingSolidBodyMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::interpolatingSolidBodyMotionSolver::calcScale()
{
    const pointMesh& pMesh = pointMesh::New(mesh());

    pointPatchDist pDist(pMesh, patchSet_, points0());

    // Scaling: 1 up to di then linear down to 0 at do away from patches
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

    // Convert the scale function to a cosine
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interpolatingSolidBodyMotionSolver::interpolatingSolidBodyMotionSolver
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    points0MotionSolver(name, mesh, dict, typeName),
    SBMFPtr_(solidBodyMotionFunction::New(coeffDict(), mesh.time())),
    patches_(wordReList(coeffDict().lookup("patches"))),
    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
    CofG_(coeffDict().lookup("CofG")),
    di_(coeffDict().lookup<scalar>("innerDistance")),
    do_(coeffDict().lookup<scalar>("outerDistance")),
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

Foam::interpolatingSolidBodyMotionSolver::~interpolatingSolidBodyMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::interpolatingSolidBodyMotionSolver::curPoints() const
{
    const pointField& points0 = this->points0();
    const septernion s = SBMFPtr_().transformation();

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
                    CofG_ + s.transformPoint(points0[pointi] - CofG_);
            }
            // Slerp septernion interpolation
            else
            {
                const septernion ss(slerp(septernion::I, s, scale_[pointi]));

                points[pointi] =
                    CofG_ + ss.transformPoint(points0[pointi] - CofG_);
            }
        }
    }

    return tpoints;
}


void Foam::interpolatingSolidBodyMotionSolver::topoChange
(
    const polyTopoChangeMap& map
)
{
    // Pending implementation of the inverse transformation of points0
    NotImplemented;
}


void Foam::interpolatingSolidBodyMotionSolver::mapMesh(const polyMeshMap& map)
{
    points0MotionSolver::mapMesh(map);

    // scale is resized by the meshToMesh mapper
    scale_ = Zero;
    calcScale();
}


// ************************************************************************* //
