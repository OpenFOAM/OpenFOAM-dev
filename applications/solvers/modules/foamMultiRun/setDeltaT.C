/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "setDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::setDeltaT(Time& runTime, const PtrList<solver>& solvers)
{
    if
    (
        runTime.timeIndex() == 0
     && runTime.controlDict().lookupOrDefault("adjustTimeStep", false)
    )
    {
        scalar deltaT = great;

        forAll(solvers, i)
        {
            if (solvers[i].transient())
            {
                deltaT = min(deltaT, solvers[i].maxDeltaT());
            }
        }

        if (deltaT != great)
        {
            runTime.setDeltaT(min(runTime.deltaTValue(), deltaT));
        }
    }
}


void Foam::adjustDeltaT(Time& runTime, const PtrList<solver>& solvers)
{
    // Update the time-step limited by the solvers maxDeltaT
    if (runTime.controlDict().lookupOrDefault("adjustTimeStep", false))
    {
        scalar deltaT = 1.2*runTime.deltaTValue();

        bool transient = false;

        forAll(solvers, i)
        {
            if (solvers[i].transient())
            {
                transient = true;
                deltaT = min(deltaT, solvers[i].maxDeltaT());
            }
        }

        if (transient)
        {
            runTime.setDeltaT(deltaT);
            Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
        }
    }
}


// ************************************************************************* //
