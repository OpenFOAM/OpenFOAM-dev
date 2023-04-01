// Calculate absolute flux from the mapped surface velocity
phic = mesh.Sf() & Ucf();

correctUphiBCs(Uc, phic, true);

fv::correctPhi
(
    phic,
    Uc,
    p,
    autoPtr<volScalarField>(),
    autoPtr<volScalarField>(),
    pressureReference,
    pimple
);

#include "continuityErrs.H"

// Make the flux relative to the mesh motion
fvc::makeRelative(phic, Uc);
