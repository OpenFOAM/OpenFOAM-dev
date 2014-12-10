{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir
    (
        IOobject
        (
            "phir",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mixture.cAlpha()*mag(phi/mesh.magSf())*mixture.nHatf()
    );

    for (int gCorr=0; gCorr<nAlphaCorr; gCorr++)
    {
        // Create the limiter to be used for all phase-fractions
        scalarField allLambda(mesh.nFaces(), 1.0);

        // Split the limiter into a surfaceScalarField
        slicedSurfaceScalarField lambda
        (
            IOobject
            (
                "lambda",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimless,
            allLambda,
            false   // Use slices for the couples
        );


        // Create the complete convection flux for alpha1
        surfaceScalarField phiAlpha1
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha3, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        // Create the bounded (upwind) flux for alpha1
        surfaceScalarField phiAlpha1BD
        (
            upwind<scalar>(mesh, phi).flux(alpha1)
        );

        // Calculate the flux correction for alpha1
        phiAlpha1 -= phiAlpha1BD;

        // Calculate the limiter for alpha1
        MULES::limiter
        (
            allLambda,
            1.0/runTime.deltaT().value(),
            geometricOneField(),
            alpha1,
            phiAlpha1BD,
            phiAlpha1,
            zeroField(),
            zeroField(),
            1,
            0,
            3
        );

        // Create the complete flux for alpha2
        surfaceScalarField phiAlpha2
        (
            fvc::flux
            (
                phi,
                alpha2,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(phir, alpha1, alpharScheme),
                alpha2,
                alpharScheme
            )
        );

        // Create the bounded (upwind) flux for alpha2
        surfaceScalarField phiAlpha2BD
        (
            upwind<scalar>(mesh, phi).flux(alpha2)
        );

        // Calculate the flux correction for alpha2
        phiAlpha2 -= phiAlpha2BD;

        // Further limit the limiter for alpha2
        MULES::limiter
        (
            allLambda,
            1.0/runTime.deltaT().value(),
            geometricOneField(),
            alpha2,
            phiAlpha2BD,
            phiAlpha2,
            zeroField(),
            zeroField(),
            1,
            0,
            3
        );

        // Construct the limited fluxes
        phiAlpha1 = phiAlpha1BD + lambda*phiAlpha1;
        phiAlpha2 = phiAlpha2BD + lambda*phiAlpha2;

        // Solve for alpha1
        solve(fvm::ddt(alpha1) + fvc::div(phiAlpha1));

        // Create the diffusion coefficients for alpha2<->alpha3
        volScalarField Dc23(D23*max(alpha3, scalar(0))*pos(alpha2));
        volScalarField Dc32(D23*max(alpha2, scalar(0))*pos(alpha3));

        // Add the diffusive flux for alpha3->alpha2
        phiAlpha2 -= fvc::interpolate(Dc32)*mesh.magSf()*fvc::snGrad(alpha1);

        // Solve for alpha2
        fvScalarMatrix alpha2Eqn
        (
            fvm::ddt(alpha2)
          + fvc::div(phiAlpha2)
          - fvm::laplacian(Dc23 + Dc32, alpha2)
        );
        alpha2Eqn.solve();

        // Construct the complete mass flux
        rhoPhi =
              phiAlpha1*(rho1 - rho3)
            + (phiAlpha2 + alpha2Eqn.flux())*(rho2 - rho3)
            + phi*rho3;

        alpha3 = 1.0 - alpha1 - alpha2;
    }

    Info<< "Air phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;

    Info<< "Liquid phase volume fraction = "
        << alpha2.weightedAverage(mesh.V()).value()
        << "  Min(alpha2) = " << min(alpha2).value()
        << "  Max(alpha2) = " << max(alpha2).value()
        << endl;
}
