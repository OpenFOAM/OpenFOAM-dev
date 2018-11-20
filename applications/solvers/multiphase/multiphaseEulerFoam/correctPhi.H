nonOrthogonalSolutionControl pcorrControl(mesh, pimple.algorithmName());

CorrectPhi
(
    U,
    phi,
    p_rgh,
    dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1),
    geometricZeroField(),
    pcorrControl
);

#include "continuityErrs.H"
