/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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

#include "componentVelocityMotionSolver.H"
#include "polyTopoChangeMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(componentVelocityMotionSolver, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::direction Foam::componentVelocityMotionSolver::cmpt
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

Foam::componentVelocityMotionSolver::componentVelocityMotionSolver
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict,
    const word& type
)
:
    motionSolver(name, mesh, dict, type),
    cmptName_(coeffDict().lookup("component")),
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

Foam::componentVelocityMotionSolver::~componentVelocityMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::componentVelocityMotionSolver::movePoints(const pointField& p)
{
    // No local data to adapt
}


void Foam::componentVelocityMotionSolver::topoChange
(
    const polyTopoChangeMap& map
)
{
    // fvMesh updates pointFields.
}


void Foam::componentVelocityMotionSolver::mapMesh
(
    const polyMeshMap& map
)
{
    pointMotionU_ == Zero;
    pointMotionU_.correctBoundaryConditions();
}


void Foam::componentVelocityMotionSolver::distribute
(
    const polyDistributionMap& map
)
{}


// ************************************************************************* //
