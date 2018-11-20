nonOrthogonalSolutionControl pcorrControl(mesh, pimple.algorithmName());

CorrectPhi
(
    Uc,
    phic,
    p,
    dimensionedScalar("rAUf", dimTime, 1),
    geometricZeroField(),
    pcorrControl
);

#include "continuityErrs.H"
