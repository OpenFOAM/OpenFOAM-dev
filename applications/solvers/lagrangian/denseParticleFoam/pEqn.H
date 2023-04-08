{
    const volScalarField rAUc(1.0/UcEqn.A());
    const volScalarField r1ADUc(1/(1 + rAUc*Dc));

    const surfaceScalarField rAUcf(fvc::interpolate(rAUc));
    const surfaceScalarField r1ADUcf(1/(1 + rAUcf*Dcf));
    const surfaceScalarField rADUcf("Dp", r1ADUcf*rAUcf);

    volVectorField HbyA(constrainHbyA(rAUc*UcEqn.H(), Uc, p));

    surfaceScalarField phiHbyAD
    (
        "phiHbyAD",
        (
            r1ADUcf
           *(
                fvc::flux(HbyA)
              + alphacf*rAUcf*fvc::ddtCorr(Uc, phic, Ucf)
            )
        )
    );

    if (p.needReference())
    {
        fvc::makeRelative(phiHbyAD, Uc);
        adjustPhi(phiHbyAD, Uc, p);
        fvc::makeAbsolute(phiHbyAD, Uc);
    }

    phiHbyAD += rADUcf*(Fgf + Dcf*phid);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, Uc, phiHbyAD, rADUcf);

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(alphacf*rADUcf, p)
         ==
            fvc::ddt(alphac)
          + fvc::div(alphacf*phiHbyAD)
        );

        pEqn.setReference
        (
            pressureReference.refCell(),
            pressureReference.refValue()
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phic = phiHbyAD - pEqn.flux()/alphacf;

            // Explicitly relax pressure for momentum corrector
            p.relax();

            Uc =
                r1ADUc
               *(
                    HbyA
                  + rAUc
                   *(
                        fvc::reconstruct
                        (
                            Fgf - pEqn.flux()/alphacf/rADUcf
                          - Dcf*(phic - phid)
                        )
                      + Dc*fvc::reconstruct(phic - phid)
                      + Fd
                     )
                );

            Uc.correctBoundaryConditions();
            fvConstraints.constrain(Uc);

            // Correct Ucf if the mesh is moving
            fvc::correctUf(Ucf, Uc, phic);

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phic, Uc);
        }
    }
}

#include "continuityErrs.H"
