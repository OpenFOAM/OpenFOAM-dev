nonOrthogonalSolutionControl pcorrControl(mesh, pimple.algorithmName());

CorrectPhi
(
    U,
    phi,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    divU,
    pcorrControl
);

#include "continuityErrs.H"
