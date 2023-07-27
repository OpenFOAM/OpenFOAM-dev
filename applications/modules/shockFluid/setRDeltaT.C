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

#include "shockFluid.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::shockFluid::setRDeltaT(const surfaceScalarField& amaxSf)
{
    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    const scalar maxCo
    (
        pimpleDict.lookupOrDefault<scalar>("maxCo", 0.8)
    );

    const scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>
        (
            "rDeltaTSmoothingCoeff",
            0.02
        )
    );

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = fvc::surfaceSum(amaxSf)()()/((2*maxCo)*mesh.V());

    // Clip to user-defined maximum and minimum time-steps
    scalar minRDeltaT = gMin(rDeltaT.primitiveField());
    if (pimpleDict.found("maxDeltaT") || minRDeltaT < rootVSmall)
    {
        const scalar clipRDeltaT = 1/pimpleDict.lookup<scalar>("maxDeltaT");
        rDeltaT.max(clipRDeltaT);
        minRDeltaT = max(minRDeltaT, clipRDeltaT);
    }
    if (pimpleDict.found("minDeltaT"))
    {
        const scalar clipRDeltaT = 1/pimpleDict.lookup<scalar>("minDeltaT");
        rDeltaT.min(clipRDeltaT);
        minRDeltaT = min(minRDeltaT, clipRDeltaT);
    }

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField()) << ", " << 1/minRDeltaT << endl;

    // Update the boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

    Info<< "Smoothed flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
}


// ************************************************************************* //
