{
    volScalarField& rDeltaT = trDeltaT();

    scalar rDeltaTSmoothingCoeff
    (
        runTime.controlDict().lookupOrDefault<scalar>
        (
            "rDeltaTSmoothingCoeff",
            0.02
        )
    );

    surfaceScalarField maxPhi("maxPhi", phi);

    forAll(phases, phasei)
    {
        maxPhi = max(maxPhi, mag(phases[phasei].phi()));
    }

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.dimensionedInternalField() = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(maxPhi)().dimensionedInternalField()
       /((2*maxCo)*mesh.V())
    );

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.internalField())
        << ", " << gMax(1/rDeltaT.internalField()) << endl;
}
