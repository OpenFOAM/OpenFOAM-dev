if (mesh.changing())
{
    forAll(U.boundaryField(), patchI)
    {
        if (U.boundaryField()[patchI].fixesValue())
        {
            U.boundaryField()[patchI].initEvaluate();
        }
    }

    forAll(U.boundaryField(), patchI)
    {
        if (U.boundaryField()[patchI].fixesValue())
        {
            U.boundaryField()[patchI].evaluate();

            phi.boundaryField()[patchI] =
                U.boundaryField()[patchI]
              & mesh.Sf().boundaryField()[patchI];
        }
    }
}

{
    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", p.dimensions(), 0.0),
        pcorrTypes
    );

    surfaceScalarField rhof(fvc::interpolate(rho, "div(phi,rho)"));
    dimensionedScalar rAUf("rAUf", dimTime, 1.0);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::ddt(rho) + fvc::div(phi*rhof)
        );

        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux()/rhof;
        }
    }
}
