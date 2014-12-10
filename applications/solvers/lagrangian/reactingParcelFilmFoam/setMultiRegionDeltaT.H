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
    Reset the timestep to maintain a constant maximum Courant numbers.
    Reduction of time-step is immediate, but increase is damped to avoid
    unstable oscillations.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    if (CoNum == -GREAT)
    {
        CoNum = SMALL;
    }

    const scalar TFactorFluid = maxCo/(CoNum + SMALL);
    const scalar TFactorFilm = maxCo/(surfaceFilm.CourantNumber() + SMALL);

    const scalar dt0 = runTime.deltaTValue();

    runTime.setDeltaT
    (
        min
        (
            dt0*min(min(TFactorFluid, TFactorFilm), 1.2),
            maxDeltaT
        )
    );
}


// ************************************************************************* //
