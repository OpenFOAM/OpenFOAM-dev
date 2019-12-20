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

            Info<< "Phase 1 " << Y1[i].name()
                << ": min = " << min(Y1[i].primitiveField())
                << ", mean = " << average(Y1[i].primitiveField())
                << ", max = " << max(Y1[i].primitiveField())
                << endl;
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

            Info<< "Phase 2 " << Y2[i].name()
                << ": min = " << min(Y2[i].primitiveField())
                << ", mean = " << average(Y2[i].primitiveField())
                << ", max = " << max(Y2[i].primitiveField())
                << endl;

        }
    }

    fluid.correctSpecies();
}
