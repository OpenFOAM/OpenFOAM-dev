{
    volVectorField HbyA(constrainHbyA(rAUc*UcEqn.H(), Uc, p));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
           fvc::flux(HbyA)
         + alphacf*rAUcf*fvc::ddtCorr(Uc, phic, Ucf)
        )
    );

    if (p.needReference())
    {
        fvc::makeRelative(phiHbyA, Uc);
        adjustPhi(phiHbyA, Uc, p);
        fvc::makeAbsolute(phiHbyA, Uc);
    }

    phiHbyA += phicForces;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, Uc, phiHbyA, rAUcf);

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(alphacf*rAUcf, p)
         ==
            fvc::ddt(alphac) + fvc::div(alphacf*phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phic = phiHbyA - pEqn.flux()/alphacf;

            p.relax();

            Uc = HbyA
               + rAUc
                *fvc::reconstruct((phicForces - pEqn.flux()/alphacf)/rAUcf);
            Uc.correctBoundaryConditions();
            fvOptions.correct(Uc);

            // Correct Ucf if the mesh is moving
            fvc::correctUf(Ucf, Uc, phic);

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phic, Uc);
        }
    }
}

#include "continuityErrs.H"
