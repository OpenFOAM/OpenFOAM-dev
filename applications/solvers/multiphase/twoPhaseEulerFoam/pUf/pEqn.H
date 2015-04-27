surfaceScalarField alpha1f("alpha1f", fvc::interpolate(alpha1));
surfaceScalarField alpha2f("alpha2f", scalar(1) - alpha1f);

surfaceScalarField alphaRho1f0
(
    "alphaRho1f0",
    fvc::interpolate
    (
        max(alpha1.oldTime(), fluid.residualAlpha(phase1))
       *rho1.oldTime()
    )
);

surfaceScalarField alphaRho2f0
(
    "alphaRho2f0",
    fvc::interpolate
    (
        max(alpha2.oldTime(), fluid.residualAlpha(phase2))
       *rho2.oldTime()
    )
);

// Drag coefficient
surfaceScalarField Kdf("Kdf", fluid.Kdf());

// Virtual-mass coefficient
surfaceScalarField Vmf("Vmf", fluid.Vmf());

// Lift and wall-lubrication forces
surfaceScalarField Ff("Ff", fluid.Ff());

tmp<surfaceScalarField> snGradAlpha1(fvc::snGrad(alpha1)*mesh.magSf());
tmp<volScalarField> D(fluid.D());

// Turbulent-dispersion force for phase 1
surfaceScalarField Ftdf1
(
    IOobject::groupName("Ftdf", phase1.name()),
    fvc::interpolate(D() + phase1.turbulence().pPrime())*snGradAlpha1()
);

// Turbulent-dispersion force for phase 2
surfaceScalarField Ftdf2
(
    IOobject::groupName("Ftdf", phase2.name()),
    fvc::interpolate(D + phase2.turbulence().pPrime())*snGradAlpha1
);


while (pimple.correct())
{
    // Update continuity errors due to temperature changes
    #include "correctContErrs.H"

    surfaceScalarField rho1f(fvc::interpolate(rho1));
    surfaceScalarField rho2f(fvc::interpolate(rho2));

    // Correct flux BCs to be consistent with the velocity BCs
    phi1.boundaryField() ==
        mrfZones.relative(mesh.Sf().boundaryField() & U1.boundaryField());
    phi2.boundaryField() ==
        mrfZones.relative(mesh.Sf().boundaryField() & U2.boundaryField());

    rAU1f =
    (
        IOobject::groupName("rAUf", phase1.name()),
        1.0
       /(
           (alphaRho1f0 + Vmf)/runTime.deltaT()
         + fvc::interpolate(U1Eqn.A())
         + Kdf
        )
    );

    rAU2f =
    (
        IOobject::groupName("rAUf", phase2.name()),
        1.0
       /(
           (alphaRho2f0 + Vmf)/runTime.deltaT()
         + fvc::interpolate(U2Eqn.A())
         + Kdf
        )
    );

    surfaceScalarField rAlphaAU1f
    (
        IOobject::groupName("rAlphaAUf", phase1.name()),
        max(alpha1f, fluid.residualAlpha(phase1))*rAU1f
    );

    surfaceScalarField rAlphaAU2f
    (
        IOobject::groupName("rAlphaAUf", phase2.name()),
        max(alpha2f, fluid.residualAlpha(phase2))*rAU2f
    );


    volScalarField rho("rho", fluid.rho());
    surfaceScalarField ghSnGradRho
    (
        "ghSnGradRho",
        ghf*fvc::snGrad(rho)*mesh.magSf()
    );

    // Add the phase-1 buoyancy force
    surfaceScalarField phiF1
    (
        IOobject::groupName("phiF", phase1.name()),
        rAlphaAU1f
       *(
            ghSnGradRho
          - alpha2f*(rho1f - rho2f)*(g & mesh.Sf())
        )
    );

    // Add the phase-2 buoyancy force
    surfaceScalarField phiF2
    (
        IOobject::groupName("phiF", phase2.name()),
        rAlphaAU2f
       *(
            ghSnGradRho
          - alpha1f*(rho2f - rho1f)*(g & mesh.Sf())
        )
    );


    // Phase-1 predicted flux
    surfaceScalarField phiHbyA1
    (
        IOobject::groupName("phiHbyA", phase1.name()),
        phi1
    );

    phiHbyA1 =
        rAU1f
       *(
            (alphaRho1f0 + Vmf)
           *mrfZones.absolute(phi1.oldTime())/runTime.deltaT()
          + (fvc::interpolate(U1Eqn.H()) & mesh.Sf())
          + Vmf*ddtPhi2
          + Kdf*mrfZones.absolute(phi2)
          - Ff
          - Ftdf1
        );

    // Phase-2 predicted flux
    surfaceScalarField phiHbyA2
    (
        IOobject::groupName("phiHbyA", phase2.name()),
        phi2
    );

    phiHbyA2 =
        rAU2f
       *(
            (alphaRho2f0 + Vmf)
           *mrfZones.absolute(phi2.oldTime())/runTime.deltaT()
          + (fvc::interpolate(U2Eqn.H()) & mesh.Sf())
          + Vmf*ddtPhi1
          + Kdf*mrfZones.absolute(phi1)
          + Ff
          + Ftdf2
       );


    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        alpha1f*(phiHbyA1 - phiF1) + alpha2f*(phiHbyA2 - phiF2)
    );
    mrfZones.makeRelative(phiHbyA);

    phiHbyA1 -= phiF1;
    phiHbyA2 -= phiF2;

    surfaceScalarField rAUf
    (
        "rAUf",
        mag(alpha1f*rAlphaAU1f + alpha2f*rAlphaAU2f)
    );

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryField(),
        (
            phiHbyA.boundaryField()
          - (
                alpha1f.boundaryField()*phi1.boundaryField()
              + alpha2f.boundaryField()*phi2.boundaryField()
            )
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    tmp<fvScalarMatrix> pEqnComp1;
    tmp<fvScalarMatrix> pEqnComp2;

    if (pimple.transonic())
    {
        surfaceScalarField phid1
        (
            IOobject::groupName("phid", phase1.name()),
            fvc::interpolate(psi1)*phi1
        );
        surfaceScalarField phid2
        (
            IOobject::groupName("phid", phase2.name()),
            fvc::interpolate(psi2)*phi2
        );

        pEqnComp1 =
            (
                contErr1
              - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
            )/rho1
          + (alpha1/rho1)*correction
            (
                psi1*fvm::ddt(p_rgh)
              + fvm::div(phid1, p_rgh) - fvm::Sp(fvc::div(phid1), p_rgh)
            );
        deleteDemandDrivenData(pEqnComp1().faceFluxCorrectionPtr());
        pEqnComp1().relax();

        pEqnComp2 =
            (
                contErr2
              - fvc::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), rho2)
            )/rho2
          + (alpha2/rho2)*correction
            (
                psi2*fvm::ddt(p_rgh)
              + fvm::div(phid2, p_rgh) - fvm::Sp(fvc::div(phid2), p_rgh)
            );
        deleteDemandDrivenData(pEqnComp2().faceFluxCorrectionPtr());
        pEqnComp2().relax();
    }
    else
    {
        pEqnComp1 =
            (
                contErr1
              - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
            )/rho1
          + (alpha1*psi1/rho1)*correction(fvm::ddt(p_rgh));

        pEqnComp2 =
            (
                contErr2
              - fvc::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), rho2)
            )/rho2
          + (alpha2*psi2/rho2)*correction(fvm::ddt(p_rgh));
    }

    // Cache p prior to solve for density update
    volScalarField p_rgh_0("p_rgh_0", p_rgh);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

        solve
        (
            pEqnComp1() + pEqnComp2() + pEqnIncomp,
            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            surfaceScalarField mSfGradp("mSfGradp", pEqnIncomp.flux()/rAUf);

            phi = phiHbyA + pEqnIncomp.flux();

            surfaceScalarField phi1s
            (
                phiHbyA1
              + rAlphaAU1f*mSfGradp
              - rAU1f*Kdf*mrfZones.absolute(phi2)
            );

            surfaceScalarField phi2s
            (
                phiHbyA2
              + rAlphaAU2f*mSfGradp
              - rAU2f*Kdf*mrfZones.absolute(phi1)
            );

            surfaceScalarField phir
            (
                ((phi2s + rAU2f*Kdf*phi1s) - (phi1s + rAU1f*Kdf*phi2s))
               /(1.0 - rAU1f*rAU2f*sqr(Kdf))
            );

            phi1 = phi - alpha2f*phir;
            phi2 = phi + alpha1f*phir;

            U1 = fvc::reconstruct(mrfZones.absolute(phi1));
            U1.correctBoundaryConditions();
            fvOptions.correct(U1);

            U2 = fvc::reconstruct(mrfZones.absolute(phi2));
            U2.correctBoundaryConditions();
            fvOptions.correct(U2);

            U = fluid.U();

            fluid.dgdt() =
            (
                alpha1*(pEqnComp2 & p_rgh)
              - alpha2*(pEqnComp1 & p_rgh)
            );
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

// Update the phase kinetic energies
K1 = 0.5*magSqr(U1);
K2 = 0.5*magSqr(U2);

// Update the pressure time-derivative if required
if (thermo1.dpdt() || thermo2.dpdt())
{
    dpdt = fvc::ddt(p);
}

#include "pUf/DDtU.H"
