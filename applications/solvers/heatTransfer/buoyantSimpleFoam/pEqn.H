{
    rho = thermo.rho();
    rho.relax();

    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    tUEqn.clear();

    surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(rho*HbyA)
    );

    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    bool closedVolume = adjustPhi(phiHbyA, U, p_rgh);

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rhorAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
        p_rghEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi = phiHbyA - p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rhorAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    p = p_rgh + rho*gh;

    // For closed-volume cases adjust the pressure level
    // to obey overall mass continuity
    if (!thermo.incompressible() && closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
        p_rgh = p - rho*gh;
    }

    rho = thermo.rho();
    rho.relax();
    Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value()
        << endl;
}
