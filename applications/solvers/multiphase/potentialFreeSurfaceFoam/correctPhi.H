nonOrthogonalSolutionControl pcorrControl(mesh, pimple.algorithmName());

CorrectPhi
(
    U,
    phi,
    p_gh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU)),
    geometricZeroField(),
    pcorrControl
);

#include "continuityErrs.H"
