/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

// Maximum flow Courant number
scalar maxCo(readScalar(pimple.dict().lookup("maxCo")));

// Maximum time scale
scalar maxDeltaT(pimple.dict().lookupOrDefault<scalar>("maxDeltaT", GREAT));

// Smoothing parameter (0-1) when smoothing iterations > 0
scalar rDeltaTSmoothingCoeff
(
    pimple.dict().lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
);

// Maximum change in cell temperature per iteration (relative to previous value)
scalar alphaTemp(pimple.dict().lookupOrDefault("alphaTemp", 0.05));


// ************************************************************************* //
