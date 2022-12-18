if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    tmp<volScalarField> trSubDeltaT;

    if (LTS)
    {
        trSubDeltaT =
            fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
    }

    // Create a temporary alphaPhi1 to accumulate the sub-cycled alphaPhi1
    tmp<surfaceScalarField> talphaPhi1
    (
        surfaceScalarField::New
        (
            "alphaPhi1",
            mesh,
            dimensionedScalar(alphaPhi1.dimensions(), 0)
        )
    );

    List<volScalarField*> alphaPtrs({&alpha1, &alpha2});

    for
    (
        subCycle<volScalarField, subCycleFields> alphaSubCycle
        (
            alphaPtrs,
            nAlphaSubCycles
        );
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        talphaPhi1.ref() += (runTime.deltaT()/totalDeltaT)*alphaPhi1;
    }

    alphaPhi1 = talphaPhi1();
}
else
{
    #include "alphaEqn.H"
}

mixture.correct();
rho == mixture.rho();

// Calculate the mass-flux from the accumulated alphaPhi1
rhoPhi = (alphaPhi1*(rho1 - rho2) + phi*rho2);
