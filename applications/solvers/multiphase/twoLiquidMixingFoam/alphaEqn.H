{
    word alphaScheme("div(phi,alpha)");

    surfaceScalarField alphaPhi
    (
        phi.name() + alpha1.name(),
        fvc::flux
        (
            phi,
            alpha1,
            alphaScheme
        )
    );

    MULES::explicitSolve
    (
        geometricOneField(),
        alpha1,
        phi,
        alphaPhi,
        oneField(),
        zeroField()
    );

    rhoPhi = alphaPhi*(rho1 - rho2) + phi*rho2;

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
