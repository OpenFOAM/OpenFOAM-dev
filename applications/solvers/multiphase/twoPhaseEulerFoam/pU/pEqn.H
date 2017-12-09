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
    surfaceScalarField Df1
    (
        fvc::interpolate
        (
            rAU1*(D + phase1.turbulence().pPrime())
        )
    );

    // Phase-2 turbulent dispersion and particle-pressure flux
    surfaceScalarField Df2
    (
        fvc::interpolate
        (
            rAU2*(D + phase2.turbulence().pPrime())
        )
    );

    // Cache the net diffusivity for implicit diffusion treatment in the
    // phase-fraction equation
    if (implicitPhasePressure)
    {
        fluid.pPrimeByA() = Df1 + Df2;
    }

    // Lift and wall-lubrication forces
    volVectorField F(fluid.F());

    // Phase-fraction face-gradient
    surfaceScalarField snGradAlpha1(fvc::snGrad(alpha1)*mesh.magSf());

    // Phase-1 dispersion, lift and wall-lubrication flux
    phiF1 = Df1*snGradAlpha1 + fvc::flux(rAU1*F);

    // Phase-1 dispersion, lift and wall-lubrication flux
    phiF2 = - Df2*snGradAlpha1 - fvc::flux(rAU2*F);
}


// --- Pressure corrector loop
while (pimple.correct())
{
    // Update continuity errors due to temperature changes
    #include "correctContErrs.H"

    volScalarField rho("rho", fluid.rho());

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;

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
    surfaceScalarField phiCorrCoeff1(pos0(alphaf1Bar - 0.99));
    surfaceScalarField phiCorrCoeff2(pos0(0.01 - alphaf1Bar));

    {
        surfaceScalarField::Boundary& phiCorrCoeff1Bf =
            phiCorrCoeff1.boundaryFieldRef();

        surfaceScalarField::Boundary& phiCorrCoeff2Bf =
            phiCorrCoeff2.boundaryFieldRef();

        forAll(mesh.boundary(), patchi)
        {
            // Set ddtPhiCorr to 0 on non-coupled boundaries
            if
            (
               !mesh.boundary()[patchi].coupled()
             || isA<cyclicAMIFvPatch>(mesh.boundary()[patchi])
            )
            {
                phiCorrCoeff1Bf[patchi] = 0;
                phiCorrCoeff2Bf[patchi] = 0;
            }
        }
    }

    // Phase-1 predicted flux
    surfaceScalarField phiHbyA1
    (
        IOobject::groupName("phiHbyA", phase1.name()),
        fvc::flux(HbyA1)
      + phiCorrCoeff1*fvc::interpolate(alpha1.oldTime()*rho1.oldTime()*rAU1)
       *(
            MRF.absolute(phi1.oldTime())
          - fvc::flux(U1.oldTime())
        )/runTime.deltaT()
      - phiF1()
      - phig1
    );

    // Phase-2 predicted flux
    surfaceScalarField phiHbyA2
    (
        IOobject::groupName("phiHbyA", phase2.name()),
        fvc::flux(HbyA2)
      + phiCorrCoeff2*fvc::interpolate(alpha2.oldTime()*rho2.oldTime()*rAU2)
       *(
            MRF.absolute(phi2.oldTime())
          - fvc::flux(U2.oldTime())
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
        p_rgh.boundaryFieldRef(),
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
          + correction
            (
                (alpha1/rho1)*
                (
                    psi1*fvm::ddt(p_rgh)
                  + fvm::div(phid1, p_rgh) - fvm::Sp(fvc::div(phid1), p_rgh)
                )
            );
        pEqnComp1.ref().relax();

        pEqnComp2 =
            (
                contErr2
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

        solve
        (
            pEqnComp1() + pEqnComp2() + pEqnIncomp,
            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        );

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

            // Compressibility correction for phase-fraction equations
            fluid.dgdt() =
            (
                alpha1*(pEqnComp2 & p_rgh)
              - alpha2*(pEqnComp1 & p_rgh)
            );

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

                U = alpha1*(Us1 + D1*U2) + alpha2*(Us2 + D2*U1);
                volVectorField Ur(((1 - D2)*Us1 - (1 - D1)*Us2)/(1 - D1*D2));

                U1 = U + alpha2*Ur;
                U1.correctBoundaryConditions();
                fvOptions.correct(U1);

                U2 = U - alpha1*Ur;
                U2.correctBoundaryConditions();
                fvOptions.correct(U2);

                U = fluid.U();
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

// Update the phase kinetic energies
K1 = 0.5*magSqr(U1);
K2 = 0.5*magSqr(U2);

// Update the pressure time-derivative if required
if (thermo1.dpdt() || thermo2.dpdt())
{
    dpdt = fvc::ddt(p);
}
