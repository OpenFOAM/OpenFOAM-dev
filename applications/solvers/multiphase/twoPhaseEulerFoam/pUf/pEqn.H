surfaceScalarField alphaf1("alphaf1", fvc::interpolate(alpha1));
surfaceScalarField alphaf2("alphaf2", scalar(1) - alphaf1);

surfaceScalarField alphaRhof10
(
    "alphaRhof10",
    fvc::interpolate
    (
        max(alpha1.oldTime(), phase1.residualAlpha())
       *rho1.oldTime()
    )
);

surfaceScalarField alphaRhof20
(
    "alphaRhof20",
    fvc::interpolate
    (
        max(alpha2.oldTime(), phase2.residualAlpha())
       *rho2.oldTime()
    )
);

// Drag coefficient
surfaceScalarField Kdf("Kdf", fluid.Kdf());

// Virtual-mass coefficient
surfaceScalarField Vmf("Vmf", fluid.Vmf());

surfaceScalarField rAUf1
(
    IOobject::groupName("rAUf", phase1.name()),
    1.0
   /(
        (alphaRhof10 + Vmf)/runTime.deltaT()
      + fvc::interpolate(U1Eqn.A())
      + Kdf
   )
);

surfaceScalarField rAUf2
(
    IOobject::groupName("rAUf", phase2.name()),
    1.0
   /(
        (alphaRhof20 + Vmf)/runTime.deltaT()
      + fvc::interpolate(U2Eqn.A())
      + Kdf
   )
);


// Turbulent dispersion, particle-pressure, lift and wall-lubrication forces
tmp<surfaceScalarField> Ff1;
tmp<surfaceScalarField> Ff2;
{
    // Turbulent-dispersion diffusivity
    volScalarField D(fluid.D());

    // Phase-1 turbulent dispersion and particle-pressure diffusivity
    surfaceScalarField Df1
    (
        fvc::interpolate(D + phase1.turbulence().pPrime())
    );

    // Phase-2 turbulent dispersion and particle-pressure diffusivity
    surfaceScalarField Df2
    (
        fvc::interpolate(D + phase2.turbulence().pPrime())
    );

    // Cache the net diffusivity for implicit diffusion treatment in the
    // phase-fraction equation
    if (implicitPhasePressure)
    {
        fluid.pPrimeByA() = rAUf1*Df1 + rAUf2*Df2;
    }

    // Lift and wall-lubrication forces
    surfaceScalarField Ff(fluid.Ff());

    // Phase-fraction face-gradient
    surfaceScalarField snGradAlpha1(fvc::snGrad(alpha1)*mesh.magSf());

    // Phase-1 dispersion, lift and wall-lubrication force
    Ff1 = Df1*snGradAlpha1 + Ff;

    // Phase-2 dispersion, lift and wall-lubrication force
    Ff2 = -Df2*snGradAlpha1 - Ff;
}


while (pimple.correct())
{
    // Update continuity errors due to temperature changes
    #include "correctContErrs.H"

    volScalarField rho("rho", fluid.rho());

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;

    surfaceScalarField rhof1(fvc::interpolate(rho1));
    surfaceScalarField rhof2(fvc::interpolate(rho2));

    // Correct fixed-flux BCs to be consistent with the velocity BCs
    MRF.correctBoundaryFlux(U1, phi1);
    MRF.correctBoundaryFlux(U2, phi2);

    surfaceScalarField alpharAUf1
    (
        IOobject::groupName("alpharAUf", phase1.name()),
        max(alphaf1, phase1.residualAlpha())*rAUf1
    );

    surfaceScalarField alpharAUf2
    (
        IOobject::groupName("alpharAUf", phase2.name()),
        max(alphaf2, phase2.residualAlpha())*rAUf2
    );

    surfaceScalarField ghSnGradRho
    (
        "ghSnGradRho",
        ghf*fvc::snGrad(rho)*mesh.magSf()
    );

    // Phase-1 buoyancy flux
    surfaceScalarField phig1
    (
        IOobject::groupName("phig", phase1.name()),
        alpharAUf1
       *(
            ghSnGradRho
          - alphaf2*(rhof1 - rhof2)*(g & mesh.Sf())
        )
    );

    // Phase-2 buoyancy flux
    surfaceScalarField phig2
    (
        IOobject::groupName("phig", phase2.name()),
        alpharAUf2
       *(
            ghSnGradRho
          - alphaf1*(rhof2 - rhof1)*(g & mesh.Sf())
        )
    );


    // Phase-1 predicted flux
    surfaceScalarField phiHbyA1
    (
        IOobject::groupName("phiHbyA", phase1.name()),
        phi1
    );

    phiHbyA1 =
        rAUf1
       *(
            (alphaRhof10 + Vmf)
           *MRF.absolute(phi1.oldTime())/runTime.deltaT()
          + fvc::flux(U1Eqn.H())
          + Vmf*ddtPhi2
          + Kdf*MRF.absolute(phi2)
          - Ff1()
        );

    // Phase-2 predicted flux
    surfaceScalarField phiHbyA2
    (
        IOobject::groupName("phiHbyA", phase2.name()),
        phi2
    );

    phiHbyA2 =
        rAUf2
       *(
            (alphaRhof20 + Vmf)
           *MRF.absolute(phi2.oldTime())/runTime.deltaT()
          + fvc::flux(U2Eqn.H())
          + Vmf*ddtPhi1
          + Kdf*MRF.absolute(phi1)
          - Ff2()
       );


    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        alphaf1*(phiHbyA1 - phig1) + alphaf2*(phiHbyA2 - phig2)
    );
    MRF.makeRelative(phiHbyA);

    phiHbyA1 -= phig1;
    phiHbyA2 -= phig2;

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
              + alpharAUf1*mSfGradp
              - rAUf1*Kdf*MRF.absolute(phi2)
            );

            surfaceScalarField phi2s
            (
                phiHbyA2
              + alpharAUf2*mSfGradp
              - rAUf2*Kdf*MRF.absolute(phi1)
            );

            surfaceScalarField phir
            (
                ((phi2s + rAUf2*Kdf*phi1s) - (phi1s + rAUf1*Kdf*phi2s))
               /(1.0 - rAUf1*rAUf2*sqr(Kdf))
            );

            phi1 = phi - alphaf2*phir;
            phi2 = phi + alphaf1*phir;

            U1 = fvc::reconstruct(MRF.absolute(phi1));
            U1.correctBoundaryConditions();
            fvOptions.correct(U1);

            U2 = fvc::reconstruct(MRF.absolute(phi2));
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
