/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "fvMeshMoversInterpolator.H"
#include "volFields.H"
#include "pointFields.H"
#include "points0MotionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshMovers
{
    defineTypeNameAndDebug(interpolator, 0);
    addToRunTimeSelectionTable(fvMeshMover, interpolator, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshMovers::interpolator::interpolator(fvMesh& mesh)
:
    fvMeshMover(mesh),
    meshCoeffs_(dict()),
    pointInterpolator_(mesh, meshCoeffs_),
    displacement_(meshCoeffs_.lookup("displacement")),
    points0_
    (
        displacement_
      ? new pointVectorField(points0MotionSolver::readPoints0(mesh))
      : nullptr
    ),
    velocityMotionCorrection_(mesh, dict())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshMovers::interpolator::~interpolator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::fvMeshMovers::interpolator::update()
{
    if (displacement_)
    {
        mesh().movePoints(points0_() + pointInterpolator_.curPointField()());
    }
    else
    {
        mesh().movePoints(pointInterpolator_.curPointField());
    }

    velocityMotionCorrection_.update();

    return true;
}


void Foam::fvMeshMovers::interpolator::topoChange(const polyTopoChangeMap&)
{
    NotImplemented;
}


void Foam::fvMeshMovers::interpolator::mapMesh(const polyMeshMap& map)
{
    if (displacement_)
    {
        points0_() == mesh().points();
    }
}


void Foam::fvMeshMovers::interpolator::distribute
(
    const polyDistributionMap&
)
{
    NotImplemented;
}


// ************************************************************************* //
