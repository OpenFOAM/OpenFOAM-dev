if (!runTime.restart() &&correctPhi)
{
    correctUphiBCs(U, phi, true);

    CorrectPhi
    (
        phi,
        U,
        p_gh,
        surfaceScalarField("rAUf", fvc::interpolate(rAU)),
        geometricZeroField(),
        pimple
    );
}
