{
    autoPtr<phaseSystem::heatTransferTable>
        heatTransferPtr(fluid.heatTransfer());

    phaseSystem::heatTransferTable&
        heatTransfer = heatTransferPtr();

    {
        tmp<fvScalarMatrix> he1Eqn(phase1.heEqn());

        if (he1Eqn.valid())
        {
            he1Eqn =
            (
                he1Eqn
             ==
               *heatTransfer[phase1.name()]
              + fvOptions(alpha1, rho1, phase1.thermo().he())
            );

            he1Eqn->relax();
            fvOptions.constrain(he1Eqn());
            he1Eqn->solve();
        }
    }

    {
        tmp<fvScalarMatrix> he2Eqn(phase2.heEqn());

        if (he2Eqn.valid())
        {
            he2Eqn =
            (
                he2Eqn
             ==
               *heatTransfer[phase2.name()]
              + fvOptions(alpha2, rho2, phase2.thermo().he())
            );

            he2Eqn->relax();
            fvOptions.constrain(he2Eqn());
            he2Eqn->solve();
        }
    }
}

fluid.correctThermo();
