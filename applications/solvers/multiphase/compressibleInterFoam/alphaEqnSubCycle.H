{
    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(mixture.cAlpha()*phic, max(phic));

    volScalarField divU(fvc::div(fvc::absolute(phi, U)));

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("0", rhoPhi.dimensions(), 0)
        );

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqn.H"
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        }

        rhoPhi = rhoPhiSum;
    }
    else
    {
        #include "alphaEqn.H"
    }
}
