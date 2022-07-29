{
    autoPtr<phaseSystem::specieTransferTable>
        specieTransferPtr(fluid.specieTransfer());

    phaseSystem::specieTransferTable&
        specieTransfer(specieTransferPtr());

    fluid.correctReactions();

    forAll(fluid.multicomponentPhases(), multicomponentPhasei)
    {
        phaseModel& phase = fluid.multicomponentPhases()[multicomponentPhasei];

        UPtrList<volScalarField>& Y = phase.YActiveRef();
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        forAll(Y, i)
        {
            fvScalarMatrix YiEqn
            (
                phase.YiEqn(Y[i])
             ==
               *specieTransfer[Y[i].name()]
              + fvModels.source(alpha, rho, Y[i])
            );

            YiEqn.relax();
            fvConstraints.constrain(YiEqn);
            YiEqn.solve("Yi");
            fvConstraints.constrain(Y[i]);
        }
    }

    fluid.correctSpecies();
}
