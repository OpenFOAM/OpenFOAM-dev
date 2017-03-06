IOobject alphaPhiHeader
(
    "alphaPhi",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

const bool alphaRestart = alphaPhiHeader.headerOk();

// MULES flux from previous time-step
surfaceScalarField alphaPhi
(
    alphaPhiHeader,
    phi*fvc::interpolate(alpha1)
);

// MULES Correction
tmp<surfaceScalarField> talphaPhiCorr0;
