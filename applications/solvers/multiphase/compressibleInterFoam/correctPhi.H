nonOrthogonalSolutionControl pcorrControl(mesh, pimple.algorithmName());

CorrectPhi
(
    U,
    phi,
    p,
    dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1),
    divU,
    pcorrControl
);

#include "continuityErrs.H"
