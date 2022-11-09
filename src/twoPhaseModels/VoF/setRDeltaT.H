{
    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    const scalar maxCo
    (
        pimpleDict.lookupOrDefault<scalar>("maxCo", 0.9)
    );

    const scalar maxAlphaCo
    (
        pimpleDict.lookupOrDefault<scalar>("maxAlphaCo", 0.2)
    );

    const scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
    );

    const label nAlphaSpreadIter
    (
        pimpleDict.lookupOrDefault<label>("nAlphaSpreadIter", 1)
    );

    const scalar alphaSpreadDiff
    (
        pimpleDict.lookupOrDefault<scalar>("alphaSpreadDiff", 0.2)
    );

    const scalar alphaSpreadMax
    (
        pimpleDict.lookupOrDefault<scalar>("alphaSpreadMax", 0.99)
    );

    const scalar alphaSpreadMin
    (
        pimpleDict.lookupOrDefault<scalar>("alphaSpreadMin", 0.01)
    );

    const label nAlphaSweepIter
    (
        pimpleDict.lookupOrDefault<label>("nAlphaSweepIter", 5)
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

    if (maxAlphaCo < maxCo)
    {
        // Further limit the reciprocal time-step
        // in the vicinity of the interface

        volScalarField alpha1Bar(fvc::average(alpha1));

        rDeltaT.ref() = max
        (
            rDeltaT(),
            pos0(alpha1Bar() - alphaSpreadMin)
           *pos0(alphaSpreadMax - alpha1Bar())
           *fvc::surfaceSum(mag(phi))()()
           /((2*maxAlphaCo)*mesh.V())
        );
    }

    // Update the boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    if (nAlphaSpreadIter > 0)
    {
        fvc::spread
        (
            rDeltaT,
            alpha1,
            nAlphaSpreadIter,
            alphaSpreadDiff,
            alphaSpreadMax,
            alphaSpreadMin
        );
    }

    if (nAlphaSweepIter > 0)
    {
        fvc::sweep(rDeltaT, alpha1, nAlphaSweepIter, alphaSpreadDiff);
    }

    Info<< "Smoothed flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        pimpleDict.found("rDeltaTDampingCoeff")
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        // Damping coefficient (1-0)
        const scalar rDeltaTDampingCoeff
        (
            pimpleDict.lookup<scalar>("rDeltaTDampingCoeff")
        );

        rDeltaT = max
        (
            rDeltaT,
            (scalar(1) - rDeltaTDampingCoeff)*rDeltaT0
        );

        Info<< "Damped flow time scale min/max = "
            << gMin(1/rDeltaT.primitiveField())
            << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
    }
}
