    surfaceScalarField ddtPhi1
    (
        (phi1 - phi1.oldTime())/runTime.deltaT()
    );

    surfaceScalarField ddtPhi2
    (
        (phi2 - phi2.oldTime())/runTime.deltaT()
    );
