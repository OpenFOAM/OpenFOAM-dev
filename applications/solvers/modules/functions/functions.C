/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "functions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(functions, 0);
    addToRunTimeSelectionTable(solver, functions, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::functions::functions(fvMesh& mesh)
:
    movingMesh(mesh)
{
    // Read the solverName from the subSolver or solver entry in controlDict
    const word solverName
    (
        runTime.controlDict().found("subSolver")
      ? runTime.controlDict().lookup("subSolver")
      : runTime.controlDict().lookup("solver")
    );

    // Instantiate the selected solver
    solverPtr = (solver::New(solverName, mesh));

    // Set all registered objects to NO_WRITE
    // so only those created by the functionObjects are written
    for
    (
        objectRegistry::iterator iter = mesh.objectRegistry::begin();
        iter != mesh.objectRegistry::end();
        ++iter
    )
    {
        iter()->writeOpt() = IOobject::NO_WRITE;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::functions::~functions()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::functions::maxDeltaT() const
{
    return solverPtr->maxDeltaT();
}


void Foam::solvers::functions::prePredictor()
{}


void Foam::solvers::functions::momentumPredictor()
{}


void Foam::solvers::functions::thermophysicalPredictor()
{}


void Foam::solvers::functions::pressureCorrector()
{}


void Foam::solvers::functions::postCorrector()
{}


void Foam::solvers::functions::postSolve()
{}


// ************************************************************************* //
