nonOrthogonalSolutionControl pcorrControl(mesh, pimple.algorithmName());

CorrectPhi
(
    U,
    phi,
    p,
    rho,
    psi,
    dimensionedScalar("rAUf", dimTime, 1),
    divrhoU(),
    pcorrControl
);
