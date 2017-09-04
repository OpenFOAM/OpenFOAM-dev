for (int Ecorr=0; Ecorr<nEnergyCorrectors; Ecorr++)
{
    fluid.correctEnergyTransport();

    autoPtr<phaseSystem::heatTransferTable>
        heatTransferPtr(fluid.heatTransfer());

    phaseSystem::heatTransferTable& heatTransfer = heatTransferPtr();

    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];

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
            fvOptions.constrain(EEqn.ref());
            EEqn->solve();
            fvOptions.correct(phase.thermo().he());
        }
    }

    fluid.correctThermo();
    fluid.correct();
}


forAll(phases, phasei)
{
    phaseModel& phase = phases[phasei];

    Info<< phase.name() << " min/max T "
        << min(phase.thermo().T()).value()
        << " - "
        << max(phase.thermo().T()).value()
        << endl;
}
