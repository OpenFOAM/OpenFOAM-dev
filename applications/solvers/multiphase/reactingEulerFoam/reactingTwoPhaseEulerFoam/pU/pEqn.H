surfaceScalarField alphaf1("alphaf1", fvc::interpolate(alpha1));
surfaceScalarField alphaf2("alphaf2", scalar(1) - alphaf1);

volScalarField rAU1
(
    IOobject::groupName("rAU", phase1.name()),
    1.0
   /(
        U1Eqn.A()
      + max(phase1.residualAlpha() - alpha1, scalar(0))
       *rho1/runTime.deltaT()
    )
);
volScalarField rAU2
(
    IOobject::groupName("rAU", phase2.name()),
    1.0
   /(
        U2Eqn.A()
      + max(phase2.residualAlpha() - alpha2, scalar(0))
       *rho2/runTime.deltaT()
    )
);

surfaceScalarField alpharAUf1
(
    fvc::interpolate(max(alpha1, phase1.residualAlpha())*rAU1)
);
surfaceScalarField alpharAUf2
(
    fvc::interpolate(max(alpha2, phase2.residualAlpha())*rAU2)
);

// Turbulent diffusion, particle-pressure, lift and wall-lubrication fluxes
tmp<surfaceScalarField> phiF1;
tmp<surfaceScalarField> phiF2;
{
    // Turbulent-dispersion diffusivity
    volScalarField D(fluid.D());

    // Phase-1 turbulent dispersion and particle-pressure flux
    tmp<surfaceScalarField> DbyA1
    (
        fvc::interpolate
        (
            rAU1*(D + phase1.turbulence().pPrime())
        )
    );

    // Phase-2 turbulent dispersion and particle-pressure flux
    tmp<surfaceScalarField> DbyA2
    (
        fvc::interpolate
        (
            rAU2*(D + phase2.turbulence().pPrime())
        )
    );

    // Lift and wall-lubrication forces
    volVectorField F(fluid.F());

    // Phase-fraction face-gradient
    surfaceScalarField snGradAlpha1(fvc::snGrad(alpha1)*mesh.magSf());

    // Phase-1 dispersion, lift and wall-lubrication flux
    phiF1 =
    (
        DbyA1()*snGradAlpha1
      + (fvc::interpolate(rAU1*F) & mesh.Sf())
    );

    // Phase-2 dispersion, lift and wall-lubrication flux
    phiF2 =
    (
      - DbyA2()*snGradAlpha1
      - (fvc::interpolate(rAU2*F) & mesh.Sf())
    );

    // Cache the phase diffusivities for implicit treatment in the
    // phase-fraction equation
    if (implicitPhasePressure)
    {
        phase1.DbyA(DbyA1);
        phase2.DbyA(DbyA2);
    }
}


// --- Pressure corrector loop
while (pimple.correct())
{
    // Update continuity errors due to temperature changes
    fluid.correct();

    // Correct fixed-flux BCs to be consistent with the velocity BCs
    MRF.correctBoundaryFlux(U1, phi1);
    MRF.correctBoundaryFlux(U2, phi2);

    volVectorField HbyA1
    (
        IOobject::groupName("HbyA", phase1.name()),
        U1
    );
    HbyA1 =
        rAU1
       *(
            U1Eqn.H()
          + max(phase1.residualAlpha() - alpha1, scalar(0))
           *rho1*U1.oldTime()/runTime.deltaT()
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
         +  max(phase2.residualAlpha() - alpha2, scalar(0))
           *rho2*U2.oldTime()/runTime.deltaT()
        );

    // Mean density for buoyancy force and p_rgh -> p
    volScalarField rho("rho", fluid.rho());

    surfaceScalarField ghSnGradRho
    (
        "ghSnGradRho",
        ghf*fvc::snGrad(rho)*mesh.magSf()
    );

    surfaceScalarField phig1
    (
        alpharAUf1
       *(
           ghSnGradRho
         - alphaf2*fvc::interpolate(rho1 - rho2)*(g & mesh.Sf())
        )
    );

    surfaceScalarField phig2
    (
        alpharAUf2
       *(
           ghSnGradRho
         - alphaf1*fvc::interpolate(rho2 - rho1)*(g & mesh.Sf())
        )
    );


    // ddtPhiCorr filter -- only apply in pure(ish) phases
    surfaceScalarField alphaf1Bar(fvc::interpolate(fvc::average(alphaf1)));
    surfaceScalarField phiCorrCoeff1(pos(alphaf1Bar - 0.99));
    surfaceScalarField phiCorrCoeff2(pos(0.01 - alphaf1Bar));

    forAll(mesh.boundary(), patchi)
    {
        // Set ddtPhiCorr to 0 on non-coupled boundaries
        if
        (
            !mesh.boundary()[patchi].coupled()
         || isA<cyclicAMIFvPatch>(mesh.boundary()[patchi])
        )
        {
            phiCorrCoeff1.boundaryField()[patchi] = 0;
            phiCorrCoeff2.boundaryField()[patchi] = 0;
        }
    }

    // Phase-1 predicted flux
    surfaceScalarField phiHbyA1
    (
        IOobject::groupName("phiHbyA", phase1.name()),
        (fvc::interpolate(HbyA1) & mesh.Sf())
      + phiCorrCoeff1*fvc::interpolate(alpha1.oldTime()*rho1.oldTime()*rAU1)
       *(
            MRF.absolute(phi1.oldTime())
          - (fvc::interpolate(U1.oldTime()) & mesh.Sf())
        )/runTime.deltaT()
      - phiF1()
      - phig1
    );

    // Phase-2 predicted flux
    surfaceScalarField phiHbyA2
    (
        IOobject::groupName("phiHbyA", phase2.name()),
        (fvc::interpolate(HbyA2) & mesh.Sf())
      + phiCorrCoeff2*fvc::interpolate(alpha2.oldTime()*rho2.oldTime()*rAU2)
       *(
            MRF.absolute(phi2.oldTime())
          - (fvc::interpolate(U2.oldTime()) & mesh.Sf())
        )/runTime.deltaT()
      - phiF2()
      - phig2
    );

    // Face-drag coefficients
    surfaceScalarField rAUKd1(fvc::interpolate(rAU1*Kd));
    surfaceScalarField rAUKd2(fvc::interpolate(rAU2*Kd));

    // Construct the mean predicted flux
    // including explicit drag contributions based on absolute fluxes
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        alphaf1*(phiHbyA1 + rAUKd1*MRF.absolute(phi2))
      + alphaf2*(phiHbyA2 + rAUKd2*MRF.absolute(phi1))
    );
    MRF.makeRelative(phiHbyA);

    // Construct pressure "diffusivity"
    surfaceScalarField rAUf
    (
        "rAUf",
        mag(alphaf1*alpharAUf1 + alphaf2*alpharAUf2)
    );

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryField(),
        (
            phiHbyA.boundaryField()
          - (
                alphaf1.boundaryField()*phi1.boundaryField()
              + alphaf2.boundaryField()*phi2.boundaryField()
            )
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    tmp<fvScalarMatrix> pEqnComp1;
    tmp<fvScalarMatrix> pEqnComp2;

    // Construct the compressibility parts of the pressure equation
    if (pimple.transonic())
    {
        if (phase1.compressible())
        {
            surfaceScalarField phid1
            (
                IOobject::groupName("phid", phase1.name()),
                fvc::interpolate(psi1)*phi1
            );

            pEqnComp1 =
                (
                    phase1.continuityError()
                  - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
                )/rho1
              + (alpha1/rho1)*correction
                (
                    psi1*fvm::ddt(p_rgh)
                  + fvm::div(phid1, p_rgh) - fvm::Sp(fvc::div(phid1), p_rgh)
                );

            deleteDemandDrivenData(pEqnComp1().faceFluxCorrectionPtr());
            pEqnComp1().relax();
        }

        if (phase2.compressible())
        {
            surfaceScalarField phid2
            (
                IOobject::groupName("phid", phase2.name()),
                fvc::interpolate(psi2)*phi2
            );

            pEqnComp2 =
                (
                    phase2.continuityError()
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
    }
    else
    {
        if (phase1.compressible())
        {
            pEqnComp1 =
                (
                    phase1.continuityError()
                  - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
                )/rho1
              + (alpha1*psi1/rho1)*correction(fvm::ddt(p_rgh));
        }

        if (phase2.compressible())
        {
            pEqnComp2 =
                (
                    phase2.continuityError()
                  - fvc::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), rho2)
                )/rho2
              + (alpha2*psi2/rho2)*correction(fvm::ddt(p_rgh));
        }
    }

    if (fluid.transfersMass())
    {
        if (pEqnComp1.valid())
        {
            pEqnComp1() -= fluid.dmdt()/rho1;
        }
        else
        {
            pEqnComp1 = fvm::Su(-fluid.dmdt()/rho1, p_rgh);
        }

        if (pEqnComp2.valid())
        {
            pEqnComp2() += fluid.dmdt()/rho2;
        }
        else
        {
            pEqnComp2 = fvm::Su(fluid.dmdt()/rho2, p_rgh);
        }
    }

    // Cache p prior to solve for density update
    volScalarField p_rgh_0(p_rgh);

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

            surfaceScalarField mSfGradp("mSfGradp", pEqnIncomp.flux()/rAUf);

            // Partial-elimination phase-flux corrector
            {
                surfaceScalarField phi1s
                (
                    phiHbyA1 + alpharAUf1*mSfGradp
                );

                surfaceScalarField phi2s
                (
                    phiHbyA2 + alpharAUf2*mSfGradp
                );

                surfaceScalarField phir
                (
                    ((phi1s + rAUKd1*phi2s) - (phi2s + rAUKd2*phi1s))
                   /(1 - rAUKd1*rAUKd2)
                );

                phi1 = phi + alphaf2*phir;
                phi2 = phi - alphaf1*phir;
            }

            // Set the phase dilatation rates
            if (phase1.compressible())
            {
                phase1.divU(-pEqnComp1 & p_rgh);
            }
            if (phase2.compressible())
            {
                phase2.divU(-pEqnComp2 & p_rgh);
            }

            // Optionally relax pressure for velocity correction
            p_rgh.relax();

            mSfGradp = pEqnIncomp.flux()/rAUf;

            // Partial-elimination phase-velocity corrector
            {
                volVectorField Us1
                (
                    HbyA1
                  + fvc::reconstruct(alpharAUf1*mSfGradp - phiF1() - phig1)
                );

                volVectorField Us2
                (
                    HbyA2
                  + fvc::reconstruct(alpharAUf2*mSfGradp - phiF2() - phig2)
                );

                volScalarField D1(rAU1*Kd);
                volScalarField D2(rAU2*Kd);

                volVectorField U(alpha1*(Us1 + D1*U2) + alpha2*(Us2 + D2*U1));
                volVectorField Ur(((1 - D2)*Us1 - (1 - D1)*Us2)/(1 - D1*D2));

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
