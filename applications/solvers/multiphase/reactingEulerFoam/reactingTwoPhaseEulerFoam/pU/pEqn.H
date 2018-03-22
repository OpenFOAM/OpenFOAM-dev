const surfaceScalarField alphaf1("alphaf1", fvc::interpolate(alpha1));
const surfaceScalarField alphaf2("alphaf2", scalar(1) - alphaf1);

PtrList<volScalarField> rAUs;
rAUs.append
(
    new volScalarField
    (
        IOobject::groupName("rAU", phase1.name()),
        1.0
       /(
            U1Eqn.A()
          + byDt(max(phase1.residualAlpha() - alpha1, scalar(0))*rho1)
        )
    )
);
rAUs.append
(
    new volScalarField
    (
        IOobject::groupName("rAU", phase2.name()),
        1.0
       /(
            U2Eqn.A()
          + byDt(max(phase2.residualAlpha() - alpha2, scalar(0))*rho2)
        )
    )
);
const volScalarField& rAU1 = rAUs[0];
const volScalarField& rAU2 = rAUs[1];

const surfaceScalarField alpharAUf1
(
    fvc::interpolate(max(alpha1, phase1.residualAlpha())*rAU1)
);
const surfaceScalarField alpharAUf2
(
    fvc::interpolate(max(alpha2, phase2.residualAlpha())*rAU2)
);

// Drag coefficients
const volScalarField Kd(fluid.Kd());
const volScalarField rAUKd1(rAU1*Kd);
const volScalarField rAUKd2(rAU2*Kd);
const surfaceScalarField phiKd1(fvc::interpolate(rAUKd1));
const surfaceScalarField phiKd2(fvc::interpolate(rAUKd2));

// Explicit force fluxes
PtrList<surfaceScalarField> phiFs(fluid.phiFs(rAUs));
const surfaceScalarField& phiF1 = phiFs[0];
const surfaceScalarField& phiF2 = phiFs[1];

// --- Pressure corrector loop
while (pimple.correct())
{
    volScalarField rho("rho", fluid.rho());

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;

    // Correct fixed-flux BCs to be consistent with the velocity BCs
    MRF.correctBoundaryFlux(U1, phi1);
    MRF.correctBoundaryFlux(U2, phi2);

    // Combined buoyancy and force fluxes
    const surfaceScalarField ghSnGradRho
    (
        "ghSnGradRho",
        ghf*fvc::snGrad(rho)*mesh.magSf()
    );

    const surfaceScalarField phigF1
    (
        alpharAUf1
       *(
           ghSnGradRho
         - alphaf2*fvc::interpolate(rho1 - rho2)*(g & mesh.Sf())
        )
      + phiF1
    );

    const surfaceScalarField phigF2
    (
        alpharAUf2
       *(
           ghSnGradRho
         - alphaf1*fvc::interpolate(rho2 - rho1)*(g & mesh.Sf())
        )
      + phiF2
    );

    // Predicted velocities
    volVectorField HbyA1
    (
        IOobject::groupName("HbyA", phase1.name()),
        U1
    );
    HbyA1 =
        rAU1
       *(
            U1Eqn.H()
          + byDt(max(phase1.residualAlpha() - alpha1, scalar(0))*rho1)
           *U1.oldTime()
        );

    volVectorField HbyA2
    (
        IOobject::groupName("HbyA", phase2.name()),
        U2
    );
    HbyA2 =
        rAU2
       *(
            U2Eqn.H()
         +  byDt(max(phase2.residualAlpha() - alpha2, scalar(0))*rho2)
           *U2.oldTime()
        );

    // Correction force fluxes
    PtrList<surfaceScalarField> ddtCorrByAs(fluid.ddtCorrByAs(rAUs));

    // Predicted fluxes
    const surfaceScalarField phiHbyA1
    (
        IOobject::groupName("phiHbyA", phase1.name()),
        fvc::flux(HbyA1) - phigF1 - ddtCorrByAs[0]
    );

    const surfaceScalarField phiHbyA2
    (
        IOobject::groupName("phiHbyA", phase2.name()),
        fvc::flux(HbyA2) - phigF2 - ddtCorrByAs[1]
    );

    ddtCorrByAs.clear();

    // Drag fluxes
    PtrList<surfaceScalarField> phiKdPhis(fluid.phiKdPhis(rAUs));
    const surfaceScalarField& phiKdPhi1 = phiKdPhis[0];
    const surfaceScalarField& phiKdPhi2 = phiKdPhis[1];

    // Total predicted flux
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        alphaf1*(phiHbyA1 - phiKdPhi1) + alphaf2*(phiHbyA2 - phiKdPhi2)
    );

    MRF.makeRelative(phiHbyA);

    phiKdPhis.clear();

    // Construct pressure "diffusivity"
    const surfaceScalarField rAUf
    (
        "rAUf",
        mag(alphaf1*alpharAUf1 + alphaf2*alpharAUf2)
    );

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryFieldRef(),
        (
            phiHbyA.boundaryField()
          - (
                alphaf1.boundaryField()*phi1.boundaryField()
              + alphaf2.boundaryField()*phi2.boundaryField()
            )
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    // Construct the compressibility parts of the pressure equation
    tmp<fvScalarMatrix> pEqnComp1, pEqnComp2;
    if (phase1.compressible())
    {
        if (pimple.transonic())
        {
            const surfaceScalarField phid1
            (
                IOobject::groupName("phid", phase1.name()),
                fvc::interpolate(psi1)*phi1
            );

            pEqnComp1 =
                (
                    fvc::ddt(alpha1, rho1) + fvc::div(phase1.alphaRhoPhi())
                  - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
                )/rho1
              + correction
                (
                    (alpha1/rho1)*
                    (
                        psi1*fvm::ddt(p_rgh)
                      + fvm::div(phid1, p_rgh) - fvm::Sp(fvc::div(phid1), p_rgh)
                    )
                );

            pEqnComp1.ref().relax();
        }
        else
        {
            pEqnComp1 =
                (
                    fvc::ddt(alpha1, rho1) + fvc::div(phase1.alphaRhoPhi())
                  - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
                )/rho1
              + (alpha1*psi1/rho1)*correction(fvm::ddt(p_rgh));
        }
    }
    if (phase2.compressible())
    {
        if (pimple.transonic())
        {
            const surfaceScalarField phid2
            (
                IOobject::groupName("phid", phase2.name()),
                fvc::interpolate(psi2)*phi2
            );

            pEqnComp2 =
                (
                    fvc::ddt(alpha2, rho2) + fvc::div(phase2.alphaRhoPhi())
                  - fvc::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), rho2)
                )/rho2
              + correction
                (
                    (alpha2/rho2)*
                    (
                        psi2*fvm::ddt(p_rgh)
                      + fvm::div(phid2, p_rgh) - fvm::Sp(fvc::div(phid2), p_rgh)
                    )
                );

            pEqnComp2.ref().relax();
        }
        else
        {
            pEqnComp2 =
                (
                    fvc::ddt(alpha2, rho2) + fvc::div(phase2.alphaRhoPhi())
                  - fvc::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), rho2)
                )/rho2
              + (alpha2*psi2/rho2)*correction(fvm::ddt(p_rgh));
        }
    }

    // Add option sources
    {
        if (fvOptions.appliesToField(rho1.name()))
        {
            tmp<fvScalarMatrix> optEqn1 = fvOptions(alpha1, rho1);
            if (pEqnComp1.valid())
            {
                pEqnComp1.ref() -= (optEqn1 & rho1)/rho1;
            }
            else
            {
                pEqnComp1 = fvm::Su(- (optEqn1 & rho1)/rho1, p_rgh);
            }
        }
        if (fvOptions.appliesToField(rho2.name()))
        {
            tmp<fvScalarMatrix> optEqn2 = fvOptions(alpha2, rho2);
            if (pEqnComp2.valid())
            {
                pEqnComp2.ref() -= (optEqn2 & rho2)/rho2;
            }
            else
            {
                pEqnComp2 = fvm::Su(- (optEqn2 & rho2)/rho2, p_rgh);
            }
        }
    }

    // Add mass transfer
    {
        PtrList<volScalarField> dmdts(fluid.dmdts());
        if (dmdts.set(0))
        {
            if (pEqnComp1.valid())
            {
                pEqnComp1.ref() -= dmdts[0]/rho1;
            }
            else
            {
                pEqnComp1 = fvm::Su(- dmdts[0]/rho1, p_rgh);
            }
        }
        if (dmdts.set(1))
        {
            if (pEqnComp2.valid())
            {
                pEqnComp2.ref() -= dmdts[1]/rho2;
            }
            else
            {
                pEqnComp2 = fvm::Su(- dmdts[1]/rho2, p_rgh);
            }
        }
    }

    // Cache p prior to solve for density update
    const volScalarField p_rgh_0(p_rgh);

    // Iterate over the pressure equation to correct for non-orthogonality
    while (pimple.correctNonOrthogonal())
    {
        // Construct the transport part of the pressure equation
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

        {
            fvScalarMatrix pEqn(pEqnIncomp);

            if (pEqnComp1.valid())
            {
                pEqn += pEqnComp1();
            }

            if (pEqnComp2.valid())
            {
                pEqn += pEqnComp2();
            }

            solve
            (
                pEqn,
                mesh.solver(p_rgh.select(pimple.finalInnerIter()))
            );
        }

        // Correct fluxes and velocities on last non-orthogonal iteration
        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqnIncomp.flux();

            surfaceScalarField mSfGradp
            (
                "mSfGradp",
                pEqnIncomp.flux()/rAUf
            );

            // Partial-elimination phase-flux corrector
            {
                const surfaceScalarField phi1s
                (
                    phiHbyA1 + alpharAUf1*mSfGradp
                );

                const surfaceScalarField phi2s
                (
                    phiHbyA2 + alpharAUf2*mSfGradp
                );

                surfaceScalarField phir
                (
                    ((phi1s + phiKd1*phi2s) - (phi2s + phiKd2*phi1s))
                   /(1 - phiKd1*phiKd2)
                );

                phi1 = phi + alphaf2*phir;
                phi2 = phi - alphaf1*phir;
            }

            // Set the phase dilatation rates
            if (pEqnComp1.valid())
            {
                phase1.divU(-pEqnComp1 & p_rgh);
            }
            if (pEqnComp2.valid())
            {
                phase2.divU(-pEqnComp2 & p_rgh);
            }

            // Optionally relax pressure for velocity correction
            p_rgh.relax();

            mSfGradp = pEqnIncomp.flux()/rAUf;

            // Partial-elimination phase-velocity corrector
            {
                const volVectorField Us1
                (
                    HbyA1
                  + fvc::reconstruct(alpharAUf1*mSfGradp - phigF1)
                );

                const volVectorField Us2
                (
                    HbyA2
                  + fvc::reconstruct(alpharAUf2*mSfGradp - phigF2)
                );

                const volVectorField U
                (
                    alpha1*(Us1 + rAUKd1*U2) + alpha2*(Us2 + rAUKd2*U1)
                );

                const volVectorField Ur
                (
                    ((1 - rAUKd2)*Us1 - (1 - rAUKd1)*Us2)/(1 - rAUKd1*rAUKd2)
                );

                U1 = U + alpha2*Ur;
                U1.correctBoundaryConditions();
                fvOptions.correct(U1);

                U2 = U - alpha1*Ur;
                U2.correctBoundaryConditions();
                fvOptions.correct(U2);
            }
        }
    }

    // Update and limit the static pressure
    p = max(p_rgh + rho*gh, pMin);

    // Limit p_rgh
    p_rgh = p - rho*gh;

    // Update densities from change in p_rgh
    rho1 += psi1*(p_rgh - p_rgh_0);
    rho2 += psi2*(p_rgh - p_rgh_0);

    // Correct p_rgh for consistency with p and the updated densities
    rho = fluid.rho();
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();
}
