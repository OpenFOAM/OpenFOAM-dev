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
        tmp<volScalarField> tRho = phase.rho();
        const volScalarField& rho = tRho();
        tmp<volVectorField> tU = phase.U();
        const volVectorField& U = tU();

        fvScalarMatrix EEqn
        (
            phase.heEqn()
         ==
           *heatTransfer[phase.name()]
          + alpha*rho*(U&g)
          + fvModels.source(alpha, rho, phase.thermoRef().he())
        );

        EEqn.relax();
        fvConstraints.constrain(EEqn);
        EEqn.solve();
        fvConstraints.constrain(phase.thermoRef().he());
    }

    fluid.correctThermo();
    fluid.correctContinuityError();
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
