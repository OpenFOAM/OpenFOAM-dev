if (!mesh.steady() && !pimple.SIMPLErho())
{
    rho = thermo.rho();
}

volScalarField rAU("rAU", 1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(rho*HbyA) + rhorAUf*fvc::ddtCorr(rho, U, phi)
);

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

const bool closedVolume = adjustPhi(phiHbyA, U, p_rgh);
const bool adjustMass = closedVolume && !thermo.incompressible();

phiHbyA += phig;

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

{
    fvScalarMatrix p_rghEqnComp
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
    );

    if (pimple.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
        );

        phiHbyA -= fvc::interpolate(psi*p_rgh)*phiHbyA/fvc::interpolate(rho);

        p_rghEqnComp += fvm::div(phid, p_rgh);
    }

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    tmp<volScalarField> psip0(mesh.steady() ? tmp<volScalarField>() : psi*p);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rhorAUf, p_rgh)
        );

        fvScalarMatrix p_rghEqn(p_rghEqnComp + p_rghEqnIncomp);

        p_rghEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi = phiHbyA + p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U = HbyA
                + rAU*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rhorAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
            K = 0.5*magSqr(U);
        }
    }

    p = p_rgh + rho*gh;

    // Thermodynamic density update
    if (!mesh.steady())
    {
        thermo.correctRho(psi*p - psip0);
    }
}

// Update pressure time derivative if needed
if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}

// Solve continuity
if (!mesh.steady())
{
    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"
}
else
{
    #include "incompressible/continuityErrs.H"
}

// Pressure limiting
const bool pLimited = pressureControl.limit(p);

// For closed-volume compressible cases adjust the pressure level
// to obey overall mass continuity
if (adjustMass)
{
    p += (initialMass - fvc::domainIntegrate(thermo.rho()))
        /fvc::domainIntegrate(psi);
    p_rgh = p - rho*gh;
}

if (adjustMass || pLimited)
{
    p.correctBoundaryConditions();
}

// Density updates
if (adjustMass || pLimited || mesh.steady() || pimple.SIMPLErho())
{
    rho = thermo.rho();
}

if (mesh.steady() && !pimple.transonic())
{
    rho.relax();
}

Info<< "Min/max rho:" << min(rho).value() << ' '
    << max(rho).value() << endl;
