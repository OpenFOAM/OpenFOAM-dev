if (ign.ignited())
{
    // progress variable
    // ~~~~~~~~~~~~~~~~~
    volScalarField c("c", scalar(1) - b);

    // Unburnt gas density
    // ~~~~~~~~~~~~~~~~~~~
    volScalarField rhou(thermo.rhou());

    // Calculate flame normal etc.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    volVectorField n("n", fvc::grad(b));

    volScalarField mgb(mag(n));

    dimensionedScalar dMgb = 1.0e-3*
        (b*c*mgb)().weightedAverage(mesh.V())
       /((b*c)().weightedAverage(mesh.V()) + small)
      + dimensionedScalar(mgb.dimensions(), small);

    mgb += dMgb;

    surfaceVectorField SfHat(mesh.Sf()/mesh.magSf());
    surfaceVectorField nfVec(fvc::interpolate(n));
    nfVec += SfHat*(fvc::snGrad(b) - (SfHat & nfVec));
    nfVec /= (mag(nfVec) + dMgb);
    surfaceScalarField nf((mesh.Sf() & nfVec));
    n /= mgb;


    #include "StCorr.H"

    // Calculate turbulent flame speed flux
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    surfaceScalarField phiSt("phiSt", fvc::interpolate(rhou*StCorr*Su*Xi)*nf);

    scalar StCoNum = max
    (
        mesh.surfaceInterpolation::deltaCoeffs()
       *mag(phiSt)/(fvc::interpolate(rho)*mesh.magSf())
    ).value()*runTime.deltaTValue();

    Info<< "Max St-Courant Number = " << StCoNum << endl;

    // Create b equation
    // ~~~~~~~~~~~~~~~~~
    fvScalarMatrix bEqn
    (
        fvm::ddt(rho, b)
      + mvConvection->fvmDiv(phi, b)
      + fvm::div(phiSt, b)
      - fvm::Sp(fvc::div(phiSt), b)
      - fvm::laplacian(turbulence->alphaEff(), b)
     ==
        fvOptions(rho, b)
    );


    // Add ignition cell contribution to b-equation
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #include "ignite.H"


    // Solve for b
    // ~~~~~~~~~~~
    bEqn.relax();

    fvOptions.constrain(bEqn);

    bEqn.solve();

    fvOptions.correct(b);

    Info<< "min(b) = " << min(b).value() << endl;


    // Calculate coefficients for Gulder's flame speed correlation
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    volScalarField up(uPrimeCoef*sqrt((2.0/3.0)*turbulence->k()));
  // volScalarField up(sqrt(mag(diag(n * n) & diag(turbulence->r()))));

    volScalarField epsilon(pow(uPrimeCoef, 3)*turbulence->epsilon());

    volScalarField tauEta(sqrt(thermo.muu()/(rhou*epsilon)));

    volScalarField Reta
    (
        up
      / (
            sqrt(epsilon*tauEta)
          + dimensionedScalar(up.dimensions(), 1e-8)
        )
    );

  // volScalarField l = 0.337*k*sqrt(k)/epsilon;
  // Reta *= max((l - dimensionedScalar(dimLength, 1.5e-3))/l, 0);

    // Calculate Xi flux
    // ~~~~~~~~~~~~~~~~~
    surfaceScalarField phiXi
    (
        phiSt
      - fvc::interpolate(fvc::laplacian(turbulence->alphaEff(), b)/mgb)*nf
      + fvc::interpolate(rho)*fvc::interpolate(Su*(1.0/Xi - Xi))*nf
    );


    // Calculate mean and turbulent strain rates
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    volVectorField Ut(U + Su*Xi*n);
    volScalarField sigmat((n & n)*fvc::div(Ut) - (n & fvc::grad(Ut) & n));

    volScalarField sigmas
    (
        ((n & n)*fvc::div(U) - (n & fvc::grad(U) & n))/Xi
      + (
            (n & n)*fvc::div(Su*n)
          - (n & fvc::grad(Su*n) & n)
        )*(Xi + scalar(1))/(2*Xi)
    );


    // Calculate the unstrained laminar flame speed
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    volScalarField Su0(unstrainedLaminarFlameSpeed()());


    // Calculate the laminar flame speed in equilibrium with the applied strain
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    volScalarField SuInf(Su0*max(scalar(1) - sigmas/sigmaExt, scalar(0.01)));

    if (SuModel == "unstrained")
    {
        Su == Su0;
    }
    else if (SuModel == "equilibrium")
    {
        Su == SuInf;
    }
    else if (SuModel == "transport")
    {
        // Solve for the strained laminar flame speed
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        volScalarField Rc
        (
            (sigmas*SuInf*(Su0 - SuInf) + sqr(SuMin)*sigmaExt)
            /(sqr(Su0 - SuInf) + sqr(SuMin))
        );

        fvScalarMatrix SuEqn
        (
            fvm::ddt(rho, Su)
          + fvm::div(phi + phiXi, Su, "div(phiXi,Su)")
          - fvm::Sp(fvc::div(phiXi), Su)
          ==
          - fvm::SuSp(-rho*Rc*Su0/Su, Su)
          - fvm::SuSp(rho*(sigmas + Rc), Su)
          + fvOptions(rho, Su)
        );

        SuEqn.relax();

        fvOptions.constrain(SuEqn);

        SuEqn.solve();

        fvOptions.correct(Su);

        // Limit the maximum Su
        // ~~~~~~~~~~~~~~~~~~~~
        Su.min(SuMax);
        Su.max(SuMin);
    }
    else
    {
        FatalError
            << args.executable() << " : Unknown Su model " << SuModel
            << abort(FatalError);
    }


    // Calculate Xi according to the selected flame wrinkling model
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (XiModel == "fixed")
    {
        // Do nothing, Xi is fixed!
    }
    else if (XiModel == "algebraic")
    {
        // Simple algebraic model for Xi based on Gulders correlation
        // with a linear correction function to give a plausible profile for Xi
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Xi == scalar(1) +
            (scalar(1) + (2*XiShapeCoef)*(scalar(0.5) - b))
           *XiCoef*sqrt(up/(Su + SuMin))*Reta;
    }
    else if (XiModel == "transport")
    {
        // Calculate Xi transport coefficients based on Gulders correlation
        // and DNS data for the rate of generation
        // with a linear correction function to give a plausible profile for Xi
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        volScalarField XiEqStar
        (
            scalar(1.001) + XiCoef*sqrt(up/(Su + SuMin))*Reta
        );

        volScalarField XiEq
        (
            scalar(1.001)
          + (
                scalar(1)
              + (2*XiShapeCoef)
               *(scalar(0.5) - min(max(b, scalar(0)), scalar(1)))
            )*(XiEqStar - scalar(1.001))
        );

        volScalarField Gstar(0.28/tauEta);
        volScalarField R(Gstar*XiEqStar/(XiEqStar - scalar(1)));
        volScalarField G(R*(XiEq - scalar(1.001))/XiEq);

        // R *= (Gstar + 2*mag(dev(symm(fvc::grad(U)))))/Gstar;

        // Solve for the flame wrinkling
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fvScalarMatrix XiEqn
        (
            fvm::ddt(rho, Xi)
          + fvm::div(phi + phiXi, Xi, "div(phiXi,Xi)")
          - fvm::Sp(fvc::div(phiXi), Xi)
         ==
            rho*R
          - fvm::Sp(rho*(R - G), Xi)
          - fvm::Sp
            (
                rho*max
                (
                    sigmat - sigmas,
                    dimensionedScalar(sigmat.dimensions(), 0)
                ),
                Xi
            )
          + fvOptions(rho, Xi)
        );

        XiEqn.relax();

        fvOptions.constrain(XiEqn);

        XiEqn.solve();

        fvOptions.correct(Xi);

        // Correct boundedness of Xi
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        Xi.max(1.0);
        Info<< "max(Xi) = " << max(Xi).value() << endl;
        Info<< "max(XiEq) = " << max(XiEq).value() << endl;
    }
    else
    {
        FatalError
            << args.executable() << " : Unknown Xi model " << XiModel
            << abort(FatalError);
    }

    Info<< "Combustion progress = "
        << 100*(scalar(1) - b)().weightedAverage(mesh.V()).value() << "%"
        << endl;

    St = Xi*Su;
}
