    Info<< "Calculating field g.h\n" << endl;
    dimensionedScalar ghRef(- mag(g)*hRef);
    volScalarField gh("gh", (g & mesh.C()) - ghRef);
    surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);
