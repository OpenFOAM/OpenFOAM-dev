/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "fluidSolutionControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidSolutionControl, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidSolutionControl::fluidSolutionControl
(
    fvMesh& mesh,
    const word& algorithmName
)
:
    nonOrthogonalSolutionControl(mesh, algorithmName),
    frozenFlow_(false),
    momentumPredictor_(true),
    transonic_(false),
    consistent_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidSolutionControl::~fluidSolutionControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::fluidSolutionControl::read()
{
    if (!nonOrthogonalSolutionControl::read())
    {
        return false;
    }

    const dictionary& solutionDict = dict();

    // The solveFluid keyword is maintained here for backwards compatibility
    frozenFlow_ = !solutionDict.lookupOrDefault<bool>("solveFluid", true);
    frozenFlow_ = solutionDict.lookupOrDefault<bool>("frozenFlow", frozenFlow_);

    momentumPredictor_ =
        solutionDict.lookupOrDefault<bool>("momentumPredictor", true);
    transonic_ = solutionDict.lookupOrDefault<bool>("transonic", false);
    consistent_ = solutionDict.lookupOrDefault<bool>("consistent", false);

    return true;
}


// ************************************************************************* //
