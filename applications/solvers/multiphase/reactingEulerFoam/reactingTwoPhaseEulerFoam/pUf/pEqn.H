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
        byDt(alphaRhof10 + Vmf)
      + fvc::interpolate(U1Eqn.A())
      + Kdf
    )
);

surfaceScalarField rAUf2
(
    IOobject::groupName("rAUf", phase2.name()),
    1.0
   /(
        byDt(alphaRhof20 + Vmf)
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
        fvc::interpolate(D + phase1.pPrime())
    );

    // Phase-2 turbulent dispersion and particle-pressure diffusivity
    surfaceScalarField Df2
    (
        fvc::interpolate(D + phase2.pPrime())
    );

    // Cache the phase diffusivities for implicit treatment in the
    // phase-fraction equation
    if (implicitPhasePressure)
    {
        phase1.DbyA(rAUf1*Df1);
        phase2.DbyA(rAUf2*Df2);
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
           *byDt(MRF.absolute(phi1.oldTime()))
          + fvc::flux(U1Eqn.H())
          + Vmf*tddtPhi2()
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
           *byDt(MRF.absolute(phi2.oldTime()))
          + fvc::flux(U2Eqn.H())
          + Vmf*tddtPhi1()
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

    // Construct the compressibility parts of the pressure equation

    if (phase1.compressible())
    {
        if (pimple.transonic())
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
                    phase1.continuityError()
                  - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
                )/rho1
              + (alpha1*psi1/rho1)*correction(fvm::ddt(p_rgh));
        }
    }
    else
    {
        pEqnComp1 = fvm::Su(-(fvOptions(alpha1, rho1)&rho1)/rho1, p_rgh);
    }

    if (phase2.compressible())
    {
        if (pimple.transonic())
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
                    phase2.continuityError()
                  - fvc::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), rho2)
                )/rho2
              + (alpha2*psi2/rho2)*correction(fvm::ddt(p_rgh));
        }
    }
    else
    {
        pEqnComp2 = fvm::Su(-(fvOptions(alpha2, rho2)&rho2)/rho2, p_rgh);
    }

    if (fluid.transfersMass())
    {
        if (pEqnComp1.valid())
        {
            pEqnComp1.ref() -= fluid.dmdt(phase1)/rho1;
        }
        else
        {
            pEqnComp1 = fvm::Su(-fluid.dmdt(phase1)/rho1, p_rgh);
        }

        if (pEqnComp2.valid())
        {
            pEqnComp2.ref() -= fluid.dmdt(phase2)/rho2;
        }
        else
        {
            pEqnComp2 = fvm::Su(-fluid.dmdt(phase2)/rho2, p_rgh);
        }
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
               /(1 - rAUf1*rAUf2*sqr(Kdf))
            );

            phi1 = phi - alphaf2*phir;
            phi2 = phi + alphaf1*phir;

            U1 = fvc::reconstruct(MRF.absolute(phi1));
            U1.correctBoundaryConditions();
            fvOptions.correct(U1);

            U2 = fvc::reconstruct(MRF.absolute(phi2));
            U2.correctBoundaryConditions();
            fvOptions.correct(U2);

            // Set the phase dilatation rates
            if (pEqnComp1.valid())
            {
                phase1.divU(-pEqnComp1 & p_rgh);
            }
            if (pEqnComp2.valid())
            {
                phase2.divU(-pEqnComp2 & p_rgh);
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
