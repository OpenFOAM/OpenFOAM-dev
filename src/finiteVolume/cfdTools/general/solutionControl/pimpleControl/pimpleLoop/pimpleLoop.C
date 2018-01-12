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

#include "pimpleLoop.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pimpleLoop, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::pimpleLoop::read()
{
    nCorrPIMPLE_ =
        control_.dict().lookupOrDefault<label>("nOuterCorrectors", 1);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pimpleLoop::pimpleLoop(const solutionControl& control)
:
    control_(control),
    nCorrPIMPLE_(-1),
    corrPIMPLE_(0),
    converged_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pimpleLoop::~pimpleLoop()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::pimpleLoop::loop(correctorConvergenceControl& convergence)
{
    read();

    ++ corrPIMPLE_;

    // Handle quit conditions first
    {
        // If all corrections have been completed then end the correction loop
        if (corrPIMPLE_ > nCorrPIMPLE_)
        {
            if (convergence.hasCorrResidualControls() && nCorrPIMPLE_ > 1)
            {
                Info<< control_.algorithmName() << ": Not converged within "
                    << nCorrPIMPLE_ << " iterations" << endl;
            }

            corrPIMPLE_ = 0;

            return false;
        }

        // If converged on the last iteration then end the correction loop
        if (converged_)
        {
            Info<< control_.algorithmName() << ": Converged in "
                << corrPIMPLE_ - 1 << " iterations" << endl;

            corrPIMPLE_ = 0;
            converged_ = false;

            return false;
        }
    }

    // If we reached here, we are doing another loop
    {
        // If convergence has been reached then set the flag so that the loop
        // exits next time around
        if (!firstPIMPLEIter() && convergence.corrCriteriaSatisfied())
        {
            Info<< control_.algorithmName() << ": Converged " << nl
                << control_.algorithmSpace() << "  Doing final iteration"
                << endl;
            converged_ = true;
        }

        // Set up the next iteration by storing the index of the solution to
        // check the convergence of
        if (firstPIMPLEIter())
        {
            convergence.resetCorrSolveIndex();
        }
        else
        {
            convergence.updateCorrSolveIndex();
        }

        // Print the number of the iteration about to take place
        if (nCorrPIMPLE_ > 1)
        {
            Info<< control_.algorithmName() << ": Iteration " << corrPIMPLE_
                << endl;
        }

        return true;
    }
}


// ************************************************************************* //
