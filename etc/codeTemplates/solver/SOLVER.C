/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "NAME.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(NAME, 0);
    addToRunTimeSelectionTable(solver, NAME, fvMesh);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// Example private function
//void Foam::solvers::NAME::functionName()
//{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::NAME::NAME(fvMesh& mesh)
:
    PARENT(mesh),

    U_
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    U(U_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::NAME::~NAME()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//void Foam::solvers::NAME::moveMesh()
//{}


//void Foam::solvers::NAME::motionCorrector()
//{}


//void Foam::solvers::NAME::preSolve()
//{}


//void Foam::solvers::NAME::prePredictor()
//{}


//void Foam::solvers::NAME::momentumTransportPredictor()
//{}


//void Foam::solvers::NAME::thermophysicalTransportPredictor()
//{}


//void Foam::solvers::NAME::momentumPredictor()
//{}


//void Foam::solvers::NAME::thermophysicalPredictor()
//{}


//void Foam::solvers::NAME::pressureCorrector()
//{}


//void Foam::solvers::NAME::momentumTransportCorrector()
//{}


//void Foam::solvers::NAME::thermophysicalTransportCorrector()
//{}


//void Foam::solvers::NAME::postSolve()
//{}


// ************************************************************************* //
