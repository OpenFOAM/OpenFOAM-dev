{
    const dictionary& pimpleDict = pimple.dict();

    scalar maxCo
    (
        pimpleDict.lookupOrDefault<scalar>("maxCo", 0.9)
    );

    scalar maxAlphaCo
    (
        pimpleDict.lookupOrDefault<scalar>("maxAlphaCo", 0.2)
    );

    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
    );

    label nAlphaSpreadIter
    (
        pimpleDict.lookupOrDefault<label>("nAlphaSpreadIter", 1)
    );

    scalar alphaSpreadDiff
    (
        pimpleDict.lookupOrDefault<scalar>("alphaSpreadDiff", 0.2)
    );

    scalar alphaSpreadMax
    (
        pimpleDict.lookupOrDefault<scalar>("alphaSpreadMax", 0.99)
    );

    scalar alphaSpreadMin
    (
        pimpleDict.lookupOrDefault<scalar>("alphaSpreadMin", 0.01)
    );

    label nAlphaSweepIter
    (
        pimpleDict.lookupOrDefault<label>("nAlphaSweepIter", 5)
    );

    scalar rDeltaTDampingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    scalar maxDeltaT
    (
        pimpleDict.lookupOrDefault<scalar>("maxDeltaT", GREAT)
    );

    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.dimensionedInternalField() = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(mag(rhoPhi))().dimensionedInternalField()
       /((2*maxCo)*mesh.V()*rho.dimensionedInternalField())
    );

    if (maxAlphaCo < maxCo)
    {
        // Further limit the reciprocal time-step
        // in the vicinity of the interface

        volScalarField alpha1Bar(fvc::average(alpha1));

        rDeltaT.dimensionedInternalField() = max
        (
            rDeltaT.dimensionedInternalField(),
            pos(alpha1Bar.dimensionedInternalField() - alphaSpreadMin)
           *pos(alphaSpreadMax - alpha1Bar.dimensionedInternalField())
           *fvc::surfaceSum(mag(phi))().dimensionedInternalField()
           /((2*maxAlphaCo)*mesh.V())
        );
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.internalField())
        << ", " << gMax(1/rDeltaT.internalField()) << endl;

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
        << gMin(1/rDeltaT.internalField())
        << ", " << gMax(1/rDeltaT.internalField()) << endl;

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

        Info<< "Damped flow time scale min/max = "
            << gMin(1/rDeltaT.internalField())
            << ", " << gMax(1/rDeltaT.internalField()) << endl;
    }

    #include "alphaControls.H"

    rSubDeltaT = rDeltaT*nAlphaSubCycles;
}
