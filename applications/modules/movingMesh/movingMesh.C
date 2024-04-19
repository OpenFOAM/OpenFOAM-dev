/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

#include "movingMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(movingMesh, 0);
    addToRunTimeSelectionTable(solver, movingMesh, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::movingMesh::movingMesh(fvMesh& mesh)
:
    solver(mesh),
    maxDeltaT_
    (
        runTime.controlDict().found("maxDeltaT")
      ? runTime.controlDict().lookup<scalar>("maxDeltaT", runTime.userUnits())
      : vGreat
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::movingMesh::~movingMesh()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::movingMesh::maxDeltaT() const
{
    return maxDeltaT_;
}


void Foam::solvers::movingMesh::preSolve()
{
    // Update the mesh for topology change, mesh to mesh mapping
    mesh_.update();
}


void Foam::solvers::movingMesh::moveMesh()
{
    if (pimple.firstIter() || pimple.moveMeshOuterCorrectors())
    {
        mesh_.move();
    }
}


void Foam::solvers::movingMesh::motionCorrector()
{}


void Foam::solvers::movingMesh::prePredictor()
{}


void Foam::solvers::movingMesh::momentumPredictor()
{}


void Foam::solvers::movingMesh::thermophysicalPredictor()
{}


void Foam::solvers::movingMesh::pressureCorrector()
{}


void Foam::solvers::movingMesh::postCorrector()
{}


void Foam::solvers::movingMesh::postSolve()
{}


// ************************************************************************* //
