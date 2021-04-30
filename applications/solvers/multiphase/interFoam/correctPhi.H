// Calculate absolute flux
// from the mapped surface velocity
phi = mesh.Sf() & Uf();

correctUphiBCs(U, phi, true);

if (divU.valid())
{
    CorrectPhi
    (
        phi,
        U,
        p_rgh,
        surfaceScalarField("rAUf", fvc::interpolate(rAU())),
        divU(),
        pressureReference,
        pimple
    );
}
else
{
    CorrectPhi
    (
        phi,
        U,
        p_rgh,
        surfaceScalarField("rAUf", fvc::interpolate(rAU())),
        geometricZeroField(),
        pressureReference,
        pimple
    );
}

#include "continuityErrs.H"

// Make the flux relative to the mesh motion
fvc::makeRelative(phi, U);
