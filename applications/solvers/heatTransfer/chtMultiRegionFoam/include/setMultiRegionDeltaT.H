/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Global
    setMultiRegionDeltaT

Description
    Reset the timestep to maintain a constant maximum courant and
    diffusion Numbers. Reduction of time-step is immediate, but
    increase is damped to avoid unstable oscillations.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    if (CoNum == -GREAT)
    {
        CoNum = SMALL;
    }

    if (DiNum == -GREAT)
    {
        DiNum = SMALL;
    }

    scalar maxDeltaTFluid = maxCo/(CoNum + SMALL);
    scalar maxDeltaTSolid = maxDi/(DiNum + SMALL);

    scalar deltaTFluid =
        min
        (
            min(maxDeltaTFluid, 1.0 + 0.1*maxDeltaTFluid),
            1.2
        );

    runTime.setDeltaT
    (
        min
        (
            min(deltaTFluid, maxDeltaTSolid)*runTime.deltaT().value(),
            maxDeltaT
        )
    );

    Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
}

// ************************************************************************* //
