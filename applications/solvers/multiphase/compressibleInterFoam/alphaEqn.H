{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir(phic*mixture.nHatf());

    for (int gCorr=0; gCorr<nAlphaCorr; gCorr++)
    {
        #include "alphaSuSp.H"

        surfaceScalarField alphaPhi1
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        MULES::explicitSolve
        (
            geometricOneField(),
            alpha1,
            phi,
            alphaPhi1,
            Sp,
            Su,
            1,
            0
        );

        surfaceScalarField rho1f(fvc::interpolate(rho1));
        surfaceScalarField rho2f(fvc::interpolate(rho2));
        rhoPhi = alphaPhi1*(rho1f - rho2f) + phi*rho2f;

        alpha2 = scalar(1) - alpha1;
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Min(" << alpha2.name() << ") = " << min(alpha2).value()
        << endl;
}
