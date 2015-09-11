Info<< "Constructing momentum equations" << endl;

PtrList<fvVectorMatrix> UEqns(fluid.phases().size());

{
    autoPtr<phaseSystem::momentumTransferTable>
        momentumTransferPtr(fluid.momentumTransfer());

    phaseSystem::momentumTransferTable&
        momentumTransfer(momentumTransferPtr());

    int phasei = 0;
    forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
    {
        phaseModel& phase = iter();
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();
        volVectorField& U = phase.U();

        UEqns.set
        (
            phasei,
            new fvVectorMatrix
            (
                phase.UEqn()
             ==
               *momentumTransfer[phase.name()]
              + fvOptions(alpha, rho, U)
            )
        );

        UEqns[phasei].relax();
        fvOptions.constrain(UEqns[phasei]);

        phasei++;
    }
}
