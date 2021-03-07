MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevSigma(U)
 ==
    fvModels.source(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvConstraints.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p_gh));

    fvConstraints.constrain(U);
}
