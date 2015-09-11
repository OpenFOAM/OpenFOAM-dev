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
            fvOptions.constrain(E1Eqn());
            E1Eqn->solve();
        }
    }

    {
        tmp<fvScalarMatrix> E2eqn(phase2.heEqn());

        if (E2eqn.valid())
        {
            E2eqn =
            (
                E2eqn
             ==
               *heatTransfer[phase2.name()]
              + alpha2*rho2*(U2&g)
              + fvOptions(alpha2, rho2, phase2.thermo().he())
            );

            E2eqn->relax();
            fvOptions.constrain(E2eqn());
            E2eqn->solve();
        }
    }
}

fluid.correctThermo();

Info<< phase1.name() << " min/max T "
    << min(phase1.thermo().T()).value()
    << " - "
    << max(phase1.thermo().T()).value()
    << endl;

Info<< phase2.name() << " min/max T "
    << min(phase2.thermo().T()).value()
    << " - "
    << max(phase2.thermo().T()).value()
    << endl;
