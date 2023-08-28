/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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

#include "pimpleNoLoopControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pimpleNoLoopControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::pimpleNoLoopControl::read()
{
    if
    (
        !pisoControl::read()
     || !singleRegionConvergenceControl::read()
     || !singleRegionCorrectorConvergenceControl::read()
    )
    {
        return false;
    }

    moveMeshOuterCorrectors_ = dict().lookupOrDefault
    (
        "moveMeshOuterCorrectors",
        false
    );

    simpleRho_ =
        dict().lookupOrDefaultBackwardsCompatible<bool>
        (
            {"simpleRho", "SIMPLErho"},
            mesh().schemes().steady()
        );

    transportPredictionFirst_ =
        dict().lookupOrDefault<bool>("transportPredictionFirst", true);

    transportCorrectionFinal_ =
        dict().lookupOrDefaultBackwardsCompatible<bool>
        (
            {"transportCorrectionFinal", "turbOnFinalIterOnly"},
            true
        );

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pimpleNoLoopControl::pimpleNoLoopControl
(
    fvMesh& mesh,
    const word& algorithmName
)
:
    pisoControl(mesh, algorithmName),
    singleRegionConvergenceControl
    (
        static_cast<singleRegionSolutionControl&>(*this)
    ),
    singleRegionCorrectorConvergenceControl
    (
        static_cast<singleRegionSolutionControl&>(*this),
        "outerCorrector"
    ),
    pimpleLoopPtr_(nullptr),
    simpleRho_(false),
    transportPredictionFirst_(true),
    transportCorrectionFinal_(true)
{
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pimpleNoLoopControl::~pimpleNoLoopControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::pimpleNoLoopControl::correct()
{
    return pisoControl::correct(finalIter());
}


bool Foam::pimpleNoLoopControl::correctNonOrthogonal()
{

    return pisoControl::correctNonOrthogonal(finalIter());
}


// ************************************************************************* //
