// Calculate absolute flux from the mapped surface velocity
phic = mesh.Sf() & Ucf();

correctUphiBCs(Uc, phic, true);

CorrectPhi
(
    phic,
    Uc,
    p,
    dimensionedScalar("rAUf", dimTime, 1),
    geometricZeroField(),
    pimple
);

#include "continuityErrs.H"

// Make the flux relative to the mesh motion
fvc::makeRelative(phic, Uc);
