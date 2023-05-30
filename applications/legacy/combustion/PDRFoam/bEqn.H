tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.schemes().div("div(phi,ft_b_ha_hau)")
    )
);

volScalarField Db("Db", rho*turbulence->nuEff());

if (ign.ignited())
{
    // Calculate the unstrained laminar flame speed
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Su = unstrainedLaminarFlameSpeed()();

    const volScalarField& Xi = flameWrinkling->Xi();

    // progress variable
    // ~~~~~~~~~~~~~~~~~
    volScalarField c("c", 1.0 - b);

    // Unburnt gas density
    // ~~~~~~~~~~~~~~~~~~~
    volScalarField rhou(thermo.rhou());

    // Calculate flame normal etc.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // volVectorField n(fvc::grad(b));
    volVectorField n(fvc::reconstruct(fvc::snGrad(b)*mesh.magSf()));

    volScalarField mgb("mgb", mag(n));

    dimensionedScalar dMgb("dMgb", mgb.dimensions(), small);

    {
        volScalarField bc(b*c);

        dMgb += 1.0e-3*
            (bc*mgb)().weightedAverage(mesh.V())
           /(bc.weightedAverage(mesh.V()) + small);
    }

    mgb += dMgb;

    surfaceVectorField Sfhat(mesh.Sf()/mesh.magSf());
    surfaceVectorField nfVec(fvc::interpolate(n));
    nfVec += Sfhat*(fvc::snGrad(b) - (Sfhat & nfVec));
    nfVec /= (mag(nfVec) + dMgb);
    surfaceScalarField nf("nf", mesh.Sf() & nfVec);
    n /= mgb;

    #include "StCorr.H"

    // Calculate turbulent flame speed flux
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    surfaceScalarField phiSt("phiSt", fvc::interpolate(rhou*StCorr*St)*nf);

    #include "StCourantNo.H"

    Db = flameWrinkling->Db();

    // Create b equation
    // ~~~~~~~~~~~~~~~~~
    fvScalarMatrix bEqn
    (
        betav*fvm::ddt(rho, b)
      + mvConvection->fvmDiv(phi, b)
      + fvm::div(phiSt, b)
      - fvm::Sp(fvc::div(phiSt), b)
      - fvm::laplacian(Db, b)
     ==
        betav*fvModels.source(rho, b)
    );


    // Add ignition cell contribution to b-equation
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #include "ignite.H"

    // Solve for b
    // ~~~~~~~~~~~
    bEqn.relax();

    fvConstraints.constrain(bEqn);

    bEqn.solve();

    fvConstraints.constrain(b);

    Info<< "min(b) = " << min(b).value() << endl;

    if (thermo.containsSpecie("ft"))
    {
        volScalarField& ft = thermo.Y("ft");

        Info<< "Combustion progress = "
            << 100*(1.0 - b)().weightedAverage(mesh.V()*ft).value() << "%"
            << endl;
    }
    else
    {
        Info<< "Combustion progress = "
            << 100*(1.0 - b)().weightedAverage(mesh.V()).value() << "%"
            << endl;
    }

    // Correct the flame-wrinkling
    flameWrinkling->correct(mvConvection);

    St = Xi*Su;
}
