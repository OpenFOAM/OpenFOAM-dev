#include "alphaControls.H"

{
    // Standard face-flux compression coefficient
    surfaceScalarField phic(cAlpha*mag(phi/mesh.magSf()));

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        surfaceScalarField rhoPhiSum("rhoPhiSum", rhoPhi);

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

    rho == alpha1*rho1 + alpha2*rho2;
}
