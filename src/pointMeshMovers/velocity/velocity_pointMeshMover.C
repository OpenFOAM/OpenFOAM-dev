/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2026 OpenFOAM Foundation
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

#include "velocity_pointMeshMover.H"
#include "polyTopoChangeMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{
    defineTypeNameAndDebug(velocity, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMeshMovers::velocity::velocity
(
    const polyMesh& mesh,
    const dictionary& dict,
    const word& type
)
:
    pointMeshMover(mesh, type),
    pointMotionU_
    (
        IOobject
        (
            "pointMotionU",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMeshMovers::velocity::~velocity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointMeshMovers::velocity::movePoints(const pointField& p)
{
    // No local data that needs adapting.
}


void Foam::pointMeshMovers::velocity::topoChange(const polyTopoChangeMap& map)
{
    pointMotionU_ == Zero;
    pointMotionU_.correctBoundaryConditions();
}


void Foam::pointMeshMovers::velocity::mapMesh(const polyMeshMap& map)
{
    pointMotionU_ == Zero;
    pointMotionU_.correctBoundaryConditions();
}


void Foam::pointMeshMovers::velocity::distribute(const polyDistributionMap&)
{}


// ************************************************************************* //
