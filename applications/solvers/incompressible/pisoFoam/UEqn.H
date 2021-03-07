// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevSigma(U)
 ==
    fvModels.source(U)
);

UEqn.relax();

fvConstraints.constrain(UEqn);

if (piso.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvConstraints.constrain(U);
}
