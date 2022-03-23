// Calculate absolute flux from the mapped surface velocity
phi = mesh.Sf() & Uf();

correctUphiBCs(U, phi, true);

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
        dimensionedScalar(p.dimensions(), 0),
        pcorrTypes
    );

    surfaceScalarField rhof(fvc::interpolate(rho, "div(phi,rho)"));
    dimensionedScalar rAUf("rAUf", dimTime, 1.0);

    mesh.schemes().setFluxRequired(pcorr.name());

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

// Make the flux relative to the mesh motion
fvc::makeRelative(phi, U);
