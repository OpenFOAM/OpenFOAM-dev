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
