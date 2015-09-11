{
    autoPtr<phaseSystem::massTransferTable>
        massTransferPtr(fluid.massTransfer());

    phaseSystem::massTransferTable&
        massTransfer(massTransferPtr());

    forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
    {
        phaseModel& phase = iter();
        PtrList<volScalarField>& Y = phase.Y();
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        forAll(Y, i)
        {
            tmp<fvScalarMatrix> YiEqn(phase.YiEqn(Y[i]));

            if (YiEqn.valid())
            {
                YiEqn =
                (
                    YiEqn
                  ==
                    *massTransfer[Y[i].name()]
                   + fvOptions(alpha, rho, Y[i])
                );

                YiEqn->relax();
                YiEqn->solve(mesh.solver("Yi"));
            }
        }
    }

    fluid.massTransfer(); // updates interfacial mass flow rates
}
