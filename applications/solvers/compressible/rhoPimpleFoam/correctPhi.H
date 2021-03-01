// Calculate absolute flux
// from the mapped surface velocity
phi = mesh.Sf() & rhoUf();

correctUphiBCs(rho, U, phi, true);

CorrectPhi
(
    phi,
    p,
    rho,
    psi,
    dimensionedScalar("rAUf", dimTime, 1),
    divrhoU(),
    pimple
);

// Make the fluxes relative to the mesh-motion
fvc::makeRelative(phi, rho, U);
