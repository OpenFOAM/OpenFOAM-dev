// Calculate absolute flux
// from the mapped surface velocity
phi = mesh.Sf() & Uf();

correctUphiBCs(U, phi, true);

CorrectPhi
(
    phi,
    U,
    p_gh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU)),
    geometricZeroField(),
    pressureReference,
    pimple
);

#include "continuityErrs.H"

// Make the flux relative to the mesh motion
fvc::makeRelative(phi, U);
