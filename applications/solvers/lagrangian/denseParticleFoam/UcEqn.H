fvVectorMatrix UcEqn
(
    fvm::ddt(alphac, Uc) + fvm::div(alphaPhic, Uc)
  - fvm::Sp(fvc::ddt(alphac) + fvc::div(alphaPhic), Uc)
  + continuousPhaseTurbulence->divDevTau(Uc)
);

UcEqn.relax();

fvConstraints.constrain(UcEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UcEqn
     ==
        fvc::reconstruct
        (
            Fgf - fvc::snGrad(p)*mesh.magSf()
          - Dcf*(phic - phid)
        )
      + Dc*fvc::reconstruct(phic - phid)
      + Fd - fvm::Sp(Dc, Uc)
    );

    fvConstraints.constrain(Uc);
}
