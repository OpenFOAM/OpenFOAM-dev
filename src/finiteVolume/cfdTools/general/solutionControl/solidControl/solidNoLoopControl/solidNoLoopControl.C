/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "solidNoLoopControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidNoLoopControl, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidNoLoopControl::solidNoLoopControl
(
    fvMesh& mesh,
    const word& algorithmName,
    const pimpleLoop& loop
)
:
    nonOrthogonalSolutionControl(mesh, algorithmName),
    singleRegionConvergenceControl
    (
        static_cast<singleRegionSolutionControl&>(*this)
    ),
    singleRegionCorrectorConvergenceControl
    (
        static_cast<singleRegionSolutionControl&>(*this),
        "outerCorrector"
    ),
    loop_(loop)
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidNoLoopControl::~solidNoLoopControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::solidNoLoopControl::read()
{
    return
        nonOrthogonalSolutionControl::read()
     && readResidualControls()
     && readCorrResidualControls();
}


bool Foam::solidNoLoopControl::isFinal() const
{
    return loop_.finalPimpleIter() && finalNonOrthogonalIter();
}


// ************************************************************************* //
