for (int Ecorr=0; Ecorr<nEnergyCorrectors; Ecorr++)
{
    fluid.correctEnergyTransport();

    autoPtr<phaseSystem::heatTransferTable>
        heatTransferPtr(fluid.heatTransfer());

    phaseSystem::heatTransferTable&
        heatTransfer = heatTransferPtr();

    {
        tmp<fvScalarMatrix> E1Eqn(phase1.heEqn());

        if (E1Eqn.valid())
        {
            E1Eqn =
            (
                E1Eqn
             ==
               *heatTransfer[phase1.name()]
              + alpha1*rho1*(U1&g)
              + fvOptions(alpha1, rho1, phase1.thermo().he())
            );

            E1Eqn->relax();
            fvOptions.constrain(E1Eqn.ref());
            E1Eqn->solve();
            fvOptions.correct(phase1.thermo().he());
        }
    }

    {
        tmp<fvScalarMatrix> E2Eqn(phase2.heEqn());

        if (E2Eqn.valid())
        {
            E2Eqn =
            (
                E2Eqn
             ==
               *heatTransfer[phase2.name()]
              + alpha2*rho2*(U2&g)
              + fvOptions(alpha2, rho2, phase2.thermo().he())
            );

            E2Eqn->relax();
            fvOptions.constrain(E2Eqn.ref());
            E2Eqn->solve();
            fvOptions.correct(phase2.thermo().he());
        }
    }

    fluid.correctThermo();
    fluid.correct();
}
