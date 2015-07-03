Info<< "Constructing momentum equations" << endl;

MRF.correctBoundaryVelocity(U1);
MRF.correctBoundaryVelocity(U2);

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVol/dimTime);

volScalarField Kd(fluid.Kd());

{
    volScalarField Vm(fluid.Vm());

    const volScalarField dmdt12(posPart(fluid.dmdt()));
    const volScalarField dmdt21(negPart(fluid.dmdt()));

    {
        U1Eqn =
        (
            fvm::ddt(alpha1, rho1, U1) + fvm::div(alphaRhoPhi1, U1)
          - fvm::Sp(phase1.continuityError() + dmdt21, U1) + dmdt21*U2
          + MRF.DDt(alpha1*rho1, U1)
          + phase1.turbulence().divDevRhoReff(U1)
         ==
          - Vm
           *(
                fvm::ddt(U1)
              + fvm::div(phi1, U1)
              - fvm::Sp(fvc::div(phi1), U1)
              - phase2.DUDt()
            )
          - MRF.DDt(Vm, U1 - U2)
          + fvOptions(alpha1, rho1, U1)
        );
        U1Eqn.relax();
        U1Eqn += fvm::Sp(Kd, U1);
        fvOptions.constrain(U1Eqn);
        U1.correctBoundaryConditions();
        fvOptions.correct(U1);
    }

    {
        U2Eqn =
        (
            fvm::ddt(alpha2, rho2, U2) + fvm::div(alphaRhoPhi2, U2)
          - fvm::Sp(phase2.continuityError() - dmdt12, U2) - dmdt12*U1
          + MRF.DDt(alpha2*rho2, U2)
          + phase2.turbulence().divDevRhoReff(U2)
         ==
          - Vm
           *(
                fvm::ddt(U2)
              + fvm::div(phi2, U2)
              - fvm::Sp(fvc::div(phi2), U2)
              - phase1.DUDt()
            )
          - MRF.DDt(Vm, U2 - U1)
          + fvOptions(alpha2, rho2, U2)
        );
        U2Eqn.relax();
        U2Eqn += fvm::Sp(Kd, U2);
        fvOptions.constrain(U2Eqn);
        U2.correctBoundaryConditions();
        fvOptions.correct(U2);
    }
}
