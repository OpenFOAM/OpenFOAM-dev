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

#include "incompressibleFluid.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::incompressibleFluid::setRDeltaT()
{
    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    const scalar maxCo
    (
        pimpleDict.lookupOrDefault<scalar>("maxCo", 0.8)
    );

    const scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.02)
    );

    const scalar rDeltaTDampingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    const scalar maxDeltaT
    (
        pimpleDict.lookupOrDefault<scalar>("maxDeltaT", great)
    );

    const scalar minDeltaT
    (
        pimpleDict.lookupOrDefault<scalar>("minDeltaT", small)
    );

    const volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Set the reciprocal time-step from the local Courant number
    // and maximum and minimum time-steps
    rDeltaT.ref() = min
    (
        1/dimensionedScalar(dimTime, minDeltaT),
        max
        (
            1/dimensionedScalar(dimTime, maxDeltaT),
            fvc::surfaceSum(mag(phi))()()
           /((2*maxCo)*mesh.V())
        )
    );

    // Update the boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    Info<< "Smoothed flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1.0
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT =
            rDeltaT0
           *max(rDeltaT/rDeltaT0, scalar(1) - rDeltaTDampingCoeff);

        Info<< "Damped flow time scale min/max = "
            << gMin(1/rDeltaT.primitiveField())
            << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
    }
}


// ************************************************************************* //
