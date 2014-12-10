rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();

volScalarField rAU(1.0/UEqn().A());
volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn().H();

if (pimple.nCorrPISO() <= 1)
{
    UEqn.clear();
}

if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            (fvc::interpolate(HbyA) & mesh.Sf())
          + fvc::ddtCorr(rAU, rho, U, phi)
        )
    );

    fvOptions.makeRelative(fvc::interpolate(psi), phid);

    volScalarField Dp("Dp", rho*rAU);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(Dp, p)
          ==
            fvOptions(psi, p, rho.name())
        );

        fvOptions.constrain(pEqn);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(HbyA) & mesh.Sf())
          + fvc::ddtCorr(rAU, rho, U, phi)
        )
    );

    fvOptions.makeRelative(fvc::interpolate(rho), phiHbyA);

    volScalarField Dp("Dp", rho*rAU);

    while (pimple.correctNonOrthogonal())
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phiHbyA)
          - fvm::laplacian(Dp, p)
          ==
            fvOptions(psi, p, rho.name())
        );

        fvOptions.constrain(pEqn);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

// Recalculate density from the relaxed pressure
rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();
Info<< "rho max/min : " << max(rho).value()
    << " " << min(rho).value() << endl;

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
