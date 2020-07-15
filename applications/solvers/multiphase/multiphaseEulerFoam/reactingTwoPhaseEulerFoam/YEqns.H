{
    autoPtr<phaseSystem::specieTransferTable>
        specieTransferPtr(fluid.specieTransfer());

    phaseSystem::specieTransferTable&
        specieTransfer(specieTransferPtr());

    fluid.correctReactions();

    if (!phase1.pure())
    {
        UPtrList<volScalarField>& Y1 = phase1.YActiveRef();

        forAll(Y1, i)
        {
            fvScalarMatrix Y1iEqn
            (
                phase1.YiEqn(Y1[i])
             ==
               *specieTransfer[Y1[i].name()]
              + fvOptions(alpha1, rho1, Y1[i])
            );

            Y1iEqn.relax();
            Y1iEqn.solve("Yi");
        }
    }

    if (!phase2.pure())
    {
        UPtrList<volScalarField>& Y2 = phase2.YActiveRef();

        forAll(Y2, i)
        {
            fvScalarMatrix Y2iEqn
            (
                phase2.YiEqn(Y2[i])
             ==
               *specieTransfer[Y2[i].name()]
              + fvOptions(alpha2, rho2, Y2[i])
            );

            Y2iEqn.relax();
            Y2iEqn.solve("Yi");
        }
    }

    fluid.correctSpecies();
}
