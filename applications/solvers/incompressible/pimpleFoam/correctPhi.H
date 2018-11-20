nonOrthogonalSolutionControl pcorrControl(mesh, pimple.algorithmName());

CorrectPhi
(
    U,
    phi,
    p,
    dimensionedScalar("rAUf", dimTime, 1),
    geometricZeroField(),
    pcorrControl
);

#include "continuityErrs.H"
