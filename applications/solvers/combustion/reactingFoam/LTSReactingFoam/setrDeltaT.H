/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

{
    Info<< "Time scales min/max:" << endl;

    // Cache old reciprocal time scale field
    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Flow time scale
    {
        rDeltaT.dimensionedInternalField() =
        (
            fvc::surfaceSum(mag(phi))().dimensionedInternalField()
           /((2*maxCo)*mesh.V()*rho.dimensionedInternalField())
        );

        // Limit the largest time scale
        rDeltaT.max(1/maxDeltaT);

        Info<< "    Flow        = "
            << gMin(1/rDeltaT.internalField()) << ", "
            << gMax(1/rDeltaT.internalField()) << endl;
    }

    // Reaction source time scale
    if (alphaTemp < 1.0)
    {
        volScalarField::DimensionedInternalField rDeltaTT
        (
            mag(reaction->Sh())/(alphaTemp*rho*thermo.Cp()*T)
        );

        Info<< "    Temperature = "
            << gMin(1/(rDeltaTT.field() + VSMALL)) << ", "
            << gMax(1/(rDeltaTT.field() + VSMALL)) << endl;

        rDeltaT.dimensionedInternalField() = max
        (
            rDeltaT.dimensionedInternalField(),
            rDeltaTT
        );
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    // Spatially smooth the time scale field
    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1.0
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT = max
        (
            rDeltaT,
            (scalar(1.0) - rDeltaTDampingCoeff)*rDeltaT0
        );
    }

    Info<< "    Overall     = "
        << gMin(1/rDeltaT.internalField())
        << ", " << gMax(1/rDeltaT.internalField()) << endl;
}

// ************************************************************************* //
