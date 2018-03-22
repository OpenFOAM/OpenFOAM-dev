{
    autoPtr<phaseSystem::massTransferTable>
        massTransferPtr(fluid.massTransfer());

    phaseSystem::massTransferTable&
        massTransfer(massTransferPtr());

    forAll(fluid.multiComponentPhases(), multiComponentPhasei)
    {
        phaseModel& phase = fluid.multiComponentPhases()[multiComponentPhasei];

        UPtrList<volScalarField>& Y = phase.YActiveRef();
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        forAll(Y, i)
        {
            fvScalarMatrix YiEqn
            (
                phase.YiEqn(Y[i])
             ==
               *massTransfer[Y[i].name()]
              + fvOptions(alpha, rho, Y[i])
            );

            YiEqn.relax();
            YiEqn.solve(mesh.solver("Yi"));
        }
    }

    fluid.massTransfer(); // updates interfacial mass flow rates
}
