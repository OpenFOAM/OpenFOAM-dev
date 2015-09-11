{
    fluid.correctEnergyTransport();

    autoPtr<phaseSystem::heatTransferTable>
        heatTransferPtr(fluid.heatTransfer());

    phaseSystem::heatTransferTable& heatTransfer = heatTransferPtr();

    forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
    {
        phaseModel& phase = iter();
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();
        const volVectorField& U = phase.U();

        tmp<fvScalarMatrix> EEqn(phase.heEqn());

        if (EEqn.valid())
        {
            EEqn =
            (
                EEqn
             ==
               *heatTransfer[phase.name()]
              + alpha*rho*(U&g)
              + fvOptions(alpha, rho, phase.thermo().he())
            );

            EEqn->relax();
            fvOptions.constrain(EEqn());
            EEqn->solve();
        }
    }
}

fluid.correctThermo();

forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
{
    phaseModel& phase = iter();

    Info<< phase.name() << " min/max T "
        << min(phase.thermo().T()).value()
        << " - "
        << max(phase.thermo().T()).value()
        << endl;
}
