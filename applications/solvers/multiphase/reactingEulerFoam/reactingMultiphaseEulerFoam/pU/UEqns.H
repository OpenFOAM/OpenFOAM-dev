Info<< "Constructing momentum equations" << endl;

PtrList<fvVectorMatrix> UEqns(fluid.phases().size());

{
    autoPtr<phaseSystem::momentumTransferTable>
        momentumTransferPtr(fluid.momentumTransfer());

    phaseSystem::momentumTransferTable&
        momentumTransfer(momentumTransferPtr());

    forAll(fluid.phases(), phasei)
    {
        phaseModel& phase = fluid.phases()[phasei];

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
    }
}
