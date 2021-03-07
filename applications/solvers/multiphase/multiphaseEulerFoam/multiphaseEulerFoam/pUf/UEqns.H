Info<< "Constructing face momentum equations" << endl;

PtrList<fvVectorMatrix> UEqns(phases.size());

{
    fluid.momentumTransfer(); // !!! Update coefficients shouldn't be necessary
                              //     This should be done on demand

    autoPtr<phaseSystem::momentumTransferTable>
        momentumTransferPtr(fluid.momentumTransferf());

    phaseSystem::momentumTransferTable&
        momentumTransfer(momentumTransferPtr());

    forAll(fluid.movingPhases(), movingPhasei)
    {
        phaseModel& phase = fluid.movingPhases()[movingPhasei];

        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();
        volVectorField& U = phase.URef();

        UEqns.set
        (
            phase.index(),
            new fvVectorMatrix
            (
                phase.UfEqn()
             ==
               *momentumTransfer[phase.name()]
              + fvModels.source(alpha, rho, U)
            )
        );

        UEqns[phase.index()].relax();
        fvConstraints.constrain(UEqns[phase.index()]);
        U.correctBoundaryConditions();
        fvConstraints.constrain(U);
    }
}
