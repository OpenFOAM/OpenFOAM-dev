/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Info<< "Time scales min/max:" << endl;


{
    // Cache old time scale field
    tmp<volScalarField> trDeltaT
    (
        new volScalarField
        (
            IOobject
            (
                "rDeltaT0",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            rDeltaT
        )
    );
    const volScalarField& rDeltaT0 = trDeltaT();


    // Flow time scale
    // ~~~~~~~~~~~~~~~
    {
        rDeltaT =
            fvc::surfaceSum
            (
                mag(phi)*mesh.deltaCoeffs()/(maxCo*mesh.magSf())
            )
           /rho;

        rDeltaT.max(1.0/maxDeltaT);

        Info<< "    Flow        = "
            << gMin(1/rDeltaT.internalField()) << ", "
            << gMax(1/rDeltaT.internalField()) << endl;
    }


    // Temperature source time scale
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        scalarField tau
        (
            runTime.deltaTValue()
           *mag
            (
                parcels.hsTrans()/(mesh.V()*runTime.deltaT())
           // + fvOptions(rho, hs)
              + combustion->Sh()()
              - fvc::div(phi, 0.5*magSqr(U), "div(phi,K)")()
            )
           /rho
        );

        tau = alphaTemp*thermo.Cp()*T/(tau + ROOTVSMALL);

        Info<< "    Temperature = " << min(maxDeltaT, gMin(tau)) << ", "
            << min(maxDeltaT, gMax(tau)) << endl;

        rDeltaT.internalField() = max(rDeltaT.internalField(), 1/tau);
    }


    // Limit rate of change of time scale
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // - reduce as much as required for flow, but limit source contributions
    const dimensionedScalar deltaTRamp("deltaTRamp", dimless, 1/(1 + 0.2));
    rDeltaT = max(rDeltaT, rDeltaT0*deltaTRamp);


    // Limit the largest time scale
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rDeltaT.max(1/maxDeltaT);


    // Spatially smooth the time scale field
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

    Info<< "    Overall     = " << min(1/rDeltaT).value()
        << ", " << max(1/rDeltaT).value() << nl << endl;
}


// ************************************************************************* //
