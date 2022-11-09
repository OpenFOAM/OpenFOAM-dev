{
    #include "alphaScheme.H"

    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;
    {
        tmp<fv::ddtScheme<scalar>> tddtAlpha
        (
            fv::ddtScheme<scalar>::New
            (
                mesh,
                mesh.schemes().ddt("ddt(alpha)")
            )
        );
        const fv::ddtScheme<scalar>& ddtAlpha = tddtAlpha();

        if
        (
            isType<fv::EulerDdtScheme<scalar>>(ddtAlpha)
         || isType<fv::localEulerDdtScheme<scalar>>(ddtAlpha)
        )
        {
            ocCoeff = 0;
        }
        else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha))
        {
            if (nAlphaSubCycles > 1)
            {
                FatalErrorInFunction
                    << "Sub-cycling is not supported "
                       "with the CrankNicolson ddt scheme"
                    << exit(FatalError);
            }

            if
            (
                alphaRestart
             || mesh.time().timeIndex() > mesh.time().startTimeIndex() + 1
            )
            {
                ocCoeff =
                    refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha)
                   .ocCoeff();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Only Euler and CrankNicolson ddt schemes are supported"
                << exit(FatalError);
        }
    }

    // Set the time blending factor, 1 for Euler
    scalar cnCoeff = 1.0/(1.0 + ocCoeff);

    tmp<surfaceScalarField> phiCN(phi);

    // Calculate the Crank-Nicolson off-centred volumetric flux
    if (ocCoeff > 0)
    {
        phiCN = surfaceScalarField::New
        (
            "phiCN",
            cnCoeff*phi + (1.0 - cnCoeff)*phi.oldTime()
        );
    }

    #include "alphaSuSp.H"

    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            (
                LTS
              ? fv::localEulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
              : fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
            )
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phiCN,
                upwind<scalar>(mesh, phiCN)
            ).fvmDiv(phiCN, alpha1)
        );

        if (divU.valid())
        {
            alpha1Eqn -= Su() + fvm::Sp(Sp() + divU(), alpha1);
        }

        alpha1Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        tmp<surfaceScalarField> talphaPhi1UD(alpha1Eqn.flux());
        alphaPhi1 = talphaPhi1UD();

        if (alphaApplyPrevCorr && talphaPhi1Corr0.valid())
        {
            Info<< "Applying the previous iteration compression flux" << endl;
            MULES::correct
            (
                geometricOneField(),
                alpha1,
                alphaPhi1,
                talphaPhi1Corr0.ref(),
                oneField(),
                zeroField()
            );

            alphaPhi1 += talphaPhi1Corr0();
        }

        // Cache the upwind-flux
        talphaPhi1Corr0 = talphaPhi1UD;

        alpha2 = 1.0 - alpha1;

        mixture.correct();
    }

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        // Split operator
        tmp<surfaceScalarField> talphaPhi1Un
        (
            fvc::flux
            (
                phiCN(),
                (cnCoeff*alpha1 + (1.0 - cnCoeff)*alpha1.oldTime())(),
                compressionScheme.rewind()
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> talphaPhi1Corr(talphaPhi1Un() - alphaPhi1);
            volScalarField alpha10("alpha10", alpha1);

            if (divU.valid())
            {
                MULES::correct
                (
                    geometricOneField(),
                    alpha1,
                    talphaPhi1Un(),
                    talphaPhi1Corr.ref(),
                    Sp(),
                    (-Sp()*alpha1)(),
                    oneField(),
                    zeroField()
                );
            }
            else
            {
                MULES::correct
                (
                    geometricOneField(),
                    alpha1,
                    talphaPhi1Un(),
                    talphaPhi1Corr.ref(),
                    oneField(),
                    zeroField()
                );
            }

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                alphaPhi1 += talphaPhi1Corr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                alphaPhi1 += 0.5*talphaPhi1Corr();
            }
        }
        else
        {
            alphaPhi1 = talphaPhi1Un;

            if (divU.valid())
            {
                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha1,
                    phiCN,
                    alphaPhi1,
                    Sp(),
                    (Su() + divU()*min(alpha1(), scalar(1)))(),
                    oneField(),
                    zeroField()
                );
            }
            else
            {
                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha1,
                    phiCN,
                    alphaPhi1,
                    oneField(),
                    zeroField()
                );
            }
        }

        alpha2 = 1.0 - alpha1;

        mixture.correct();
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        talphaPhi1Corr0 = alphaPhi1 - talphaPhi1Corr0;
        talphaPhi1Corr0.ref().rename("alphaPhi1Corr0");
    }
    else
    {
        talphaPhi1Corr0.clear();
    }

    #include "rhofs.H"

    if
    (
        word(mesh.schemes().ddt("ddt(rho,U)"))
     == fv::EulerDdtScheme<vector>::typeName
     || word(mesh.schemes().ddt("ddt(rho,U)"))
     == fv::localEulerDdtScheme<vector>::typeName
    )
    {
        rhoPhi = alphaPhi1*(rho1f - rho2f) + phiCN*rho2f;
    }
    else
    {
        if (ocCoeff > 0)
        {
            // Calculate the end-of-time-step alpha flux
            alphaPhi1 =
                (alphaPhi1 - (1.0 - cnCoeff)*alphaPhi1.oldTime())/cnCoeff;
        }

        // Calculate the end-of-time-step mass flux
        rhoPhi = alphaPhi1*(rho1f - rho2f) + phi*rho2f;
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
