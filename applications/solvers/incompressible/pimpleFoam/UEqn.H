// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);

UEqn().relax();

fvOptions.constrain(UEqn());

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
}
