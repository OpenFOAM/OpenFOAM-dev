fvVectorMatrix UcEqn
(
    fvm::ddt(alphac, Uc) + fvm::div(alphaPhic, Uc)
  - fvm::Sp(fvc::ddt(alphac) + fvc::div(alphaPhic), Uc)
  + continuousPhaseTurbulence->divDevTau(Uc)
 ==
    (1.0/rhoc)*cloudSU
);

UcEqn.relax();

fvConstraints.constrain(UcEqn);

volScalarField rAUc(1.0/UcEqn.A());
volScalarField rASpUc(1.0/(UcEqn.A() - cloudSUp/rhoc));
surfaceScalarField rASpUcf("Dp", fvc::interpolate(rASpUc));

surfaceScalarField phicSUSu
(
    fvc::flux(rASpUc*cloudSUu/rhoc)
  + rASpUcf*(g & mesh.Sf())
);
surfaceScalarField phicSUSp
(
    fvc::interpolate(rASpUc*cloudSUp/rhoc)
);

if (pimple.momentumPredictor())
{
    solve
    (
        UcEqn
     ==
        fvc::reconstruct
        (
            (phicSUSu + phicSUSp*phic)/rASpUcf
          - fvc::snGrad(p)*mesh.magSf()
        )
      + (1.0/rhoc)*(fvm::Sp(cloudSUp, Uc) - cloudSUp*Uc)
    );

    fvConstraints.constrain(Uc);
}
