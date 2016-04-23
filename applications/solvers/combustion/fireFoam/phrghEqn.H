if (pimple.dict().lookupOrDefault<bool>("hydrostaticInitialization", false))
{
    volScalarField& ph_rgh = regIOobject::store
    (
        new volScalarField
        (
            IOobject
            (
                "ph_rgh",
                "0",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    );

    if (equal(runTime.value(), 0))
    {
        p = ph_rgh + rho*gh + pRef;
        thermo.correct();
        rho = thermo.rho();

        label nCorr
        (
            pimple.dict().lookupOrDefault<label>("nHydrostaticCorrectors", 5)
        );

        for (label i=0; i<nCorr; i++)
        {
            surfaceScalarField rhof("rhof", fvc::interpolate(rho));

            surfaceScalarField phig
            (
                "phig",
               -rhof*ghf*fvc::snGrad(rho)*mesh.magSf()
            );

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(ph_rgh, rho, U, phig, rhof);

            fvScalarMatrix ph_rghEqn
            (
                fvm::laplacian(rhof, ph_rgh) == fvc::div(phig)
            );

            ph_rghEqn.solve();

            p = ph_rgh + rho*gh + pRef;
            thermo.correct();
            rho = thermo.rho();

            Info<< "Hydrostatic pressure variation "
                << (max(ph_rgh) - min(ph_rgh)).value() << endl;
        }

        ph_rgh.write();

        p_rgh = ph_rgh;
    }
}
