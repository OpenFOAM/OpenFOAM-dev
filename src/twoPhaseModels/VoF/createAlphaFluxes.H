typeIOobject<surfaceScalarField> alphaPhi1Header
(
    IOobject::groupName("alphaPhi", alpha1.group()),
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

const bool alphaRestart = alphaPhi1Header.headerOk();

if (alphaRestart)
{
    Info << "Restarting alpha" << endl;
}

// Phase-1 volumetric flux
surfaceScalarField alphaPhi1
(
    alphaPhi1Header,
    phi*fvc::interpolate(alpha1)
);

// MULES Correction
tmp<surfaceScalarField> talphaPhi1Corr0;
