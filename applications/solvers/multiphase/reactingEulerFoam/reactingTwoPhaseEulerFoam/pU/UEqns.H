Info<< "Constructing momentum equations" << endl;

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVol/dimTime);

{
    autoPtr<phaseSystem::momentumTransferTable>
        momentumTransferPtr(fluid.momentumTransfer());

    phaseSystem::momentumTransferTable&
        momentumTransfer(momentumTransferPtr());

    {
        U1Eqn =
        (
            phase1.UEqn()
         ==
           *momentumTransfer[phase1.name()]
          + fvOptions(alpha1, rho1, U1)
        );
        U1Eqn.relax();
        fvOptions.constrain(U1Eqn);
        fvOptions.correct(U1);
    }

    {
        U2Eqn =
        (
            phase2.UEqn()
         ==
           *momentumTransfer[phase2.name()]
          + fvOptions(alpha2, rho2, U2)
        );
        U2Eqn.relax();
        fvOptions.constrain(U2Eqn);
        fvOptions.correct(U2);
    }
}
