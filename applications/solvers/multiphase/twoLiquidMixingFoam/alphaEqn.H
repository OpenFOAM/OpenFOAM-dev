{
    word alphaScheme("div(phi,alpha)");

    surfaceScalarField phiAlpha
    (
        fvc::flux
        (
            phi,
            alpha1,
            alphaScheme
        )
    );

    MULES::explicitSolve(alpha1, phi, phiAlpha, 1, 0);

    rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2;

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
