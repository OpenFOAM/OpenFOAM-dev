{
    autoPtr<phaseSystem::massTransferTable>
        massTransferPtr(fluid.massTransfer());

    phaseSystem::massTransferTable&
        massTransfer(massTransferPtr());

    PtrList<volScalarField>& Y1 = phase1.Y();
    PtrList<volScalarField>& Y2 = phase2.Y();

    forAll(Y1, i)
    {
        tmp<fvScalarMatrix> Y1iEqn(phase1.YiEqn(Y1[i]));

        if (Y1iEqn.valid())
        {
            Y1iEqn =
            (
                Y1iEqn
             ==
               *massTransfer[Y1[i].name()]
              + fvOptions(alpha1, rho1, Y1[i])
            );

            Y1iEqn->relax();
            Y1iEqn->solve(mesh.solver("Yi"));
        }
    }

    forAll(Y2, i)
    {
        tmp<fvScalarMatrix> Y2iEqn(phase2.YiEqn(Y2[i]));

        if (Y2iEqn.valid())
        {
            Y2iEqn =
            (
                Y2iEqn
             ==
               *massTransfer[Y2[i].name()]
              + fvOptions(alpha2, rho2, Y2[i])
            );

            Y2iEqn->relax();
            Y2iEqn->solve(mesh.solver("Yi"));
        }
    }

    fluid.massTransfer(); // updates interfacial mass flow rates
}
