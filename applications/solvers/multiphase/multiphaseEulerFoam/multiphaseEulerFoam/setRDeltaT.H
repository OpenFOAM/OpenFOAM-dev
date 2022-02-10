{
    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    const scalar maxCo
    (
        pimpleDict.lookupOrDefault<scalar>("maxCo", 0.2)
    );

    const scalar maxDeltaT
    (
        pimpleDict.lookupOrDefault<scalar>("maxDeltaT", great)
    );

    const scalar minDeltaT
    (
        pimpleDict.lookupOrDefault<scalar>("minDeltaT", small)
    );

    const scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.02)
    );

    surfaceScalarField maxPhi("maxPhi", phi);

    forAll(phases, phasei)
    {
        maxPhi = max(maxPhi, mag(phases[phasei].phi()));
    }

    // Set the reciprocal time-step from the local Courant number
    // and maximum and minimum time-steps
    rDeltaT.ref() = min
    (
        1/dimensionedScalar(dimTime, minDeltaT),
        max
        (
            1/dimensionedScalar(dimTime, maxDeltaT),
            fvc::surfaceSum(maxPhi)()()
           /((2*maxCo)*mesh.V())
        )
    );

    // Update the boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
}
