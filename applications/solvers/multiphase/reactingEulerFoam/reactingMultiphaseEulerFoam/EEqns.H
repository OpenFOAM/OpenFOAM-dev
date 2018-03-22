for (int Ecorr=0; Ecorr<nEnergyCorrectors; Ecorr++)
{
    fluid.correctEnergyTransport();

    autoPtr<phaseSystem::heatTransferTable>
        heatTransferPtr(fluid.heatTransfer());

    phaseSystem::heatTransferTable& heatTransfer = heatTransferPtr();

    forAll(fluid.anisothermalPhases(), anisothermalPhasei)
    {
        phaseModel& phase = fluid.anisothermalPhases()[anisothermalPhasei];

        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();
        const volVectorField& U = phase.U();

        fvScalarMatrix EEqn
        (
            phase.heEqn()
         ==
           *heatTransfer[phase.name()]
          + alpha*rho*(U&g)
          + fvOptions(alpha, rho, phase.thermoRef().he())
        );

        EEqn.relax();
        fvOptions.constrain(EEqn);
        EEqn.solve();
        fvOptions.correct(phase.thermoRef().he());
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
