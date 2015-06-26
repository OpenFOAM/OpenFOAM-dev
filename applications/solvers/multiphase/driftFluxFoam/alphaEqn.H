{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha1)
        );

        solve(alpha1Eqn);

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        tmp<surfaceScalarField> tphiAlphaUD(alpha1Eqn.flux());
        phiAlpha = tphiAlphaUD();

        if (alphaApplyPrevCorr && tphiAlphaCorr0.valid())
        {
            Info<< "Applying the previous iteration correction flux" << endl;

            MULES::correct
            (
                alpha1,
                phiAlpha,
                tphiAlphaCorr0(),
                mixture.alphaMax(),
                0
            );

            phiAlpha += tphiAlphaCorr0();
        }

        // Cache the upwind-flux
        tphiAlphaCorr0 = tphiAlphaUD;
    }

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        tmp<surfaceScalarField> tphiAlphaUn
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                phir,
                alpha1,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> tphiAlphaCorr(tphiAlphaUn() - phiAlpha);
            volScalarField alpha10("alpha10", alpha1);

            MULES::correct
            (
                alpha1,
                tphiAlphaUn(),
                tphiAlphaCorr(),
                mixture.alphaMax(),
                0
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                phiAlpha += tphiAlphaCorr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                phiAlpha += 0.5*tphiAlphaCorr();
            }
        }
        else
        {
            phiAlpha = tphiAlphaUn;

            MULES::explicitSolve
            (
                alpha1,
                phi,
                phiAlpha,
                mixture.alphaMax(),
                0
            );
        }
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        tphiAlphaCorr0 = phiAlpha - tphiAlphaCorr0;
    }

    alpha2 = 1.0 - alpha1;

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
