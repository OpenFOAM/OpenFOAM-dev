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

#include "velocityComponent_pointMeshMover.H"
#include "polyTopoChangeMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{
    defineTypeNameAndDebug(velocityComponent, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::direction Foam::pointMeshMovers::velocityComponent::cmpt
(
    const word& cmptName
) const
{
    if (cmptName == "x")
    {
        return vector::X;
    }
    else if (cmptName == "y")
    {
        return vector::Y;
    }
    else if (cmptName == "z")
    {
        return vector::Z;
    }
    else
    {
        FatalErrorInFunction
            << "Given component name " << cmptName << " should be x, y or z"
            << exit(FatalError);

        return 0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMeshMovers::velocityComponent::velocityComponent
(
    const polyMesh& mesh,
    const dictionary& dict,
    const word& type
)
:
    pointMeshMover(mesh, type),
    cmptName_(dict.lookup("component")),
    cmpt_(cmpt(cmptName_)),
    pointMotionU_
    (
        IOobject
        (
            "pointMotionU" + cmptName_,
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMeshMovers::velocityComponent::~velocityComponent()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointMeshMovers::velocityComponent::movePoints(const pointField& p)
{
    // No local data to adapt
}


void Foam::pointMeshMovers::velocityComponent::topoChange
(
    const polyTopoChangeMap& map
)
{
    pointMotionU_ == Zero;
    pointMotionU_.correctBoundaryConditions();
}


void Foam::pointMeshMovers::velocityComponent::mapMesh
(
    const polyMeshMap& map
)
{
    pointMotionU_ == Zero;
    pointMotionU_.correctBoundaryConditions();
}


void Foam::pointMeshMovers::velocityComponent::distribute
(
    const polyDistributionMap& map
)
{}


// ************************************************************************* //
