Info<< "Constructing face momentum equations" << endl;

PtrList<fvVectorMatrix> UEqns(phases.size());

{
    fluid.momentumTransfer(); // !!! Update coefficients shouldn't be necessary
                              //     This should be done on demand

    autoPtr<phaseSystem::momentumTransferTable>
        momentumTransferPtr(fluid.momentumTransferf());

    phaseSystem::momentumTransferTable&
        momentumTransfer(momentumTransferPtr());

    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];

        const volScalarField& alpha = phase;
        volScalarField& rho = phase.thermo().rho();
        volVectorField& U = phase.U();

        UEqns.set
        (
            phasei,
            new fvVectorMatrix
            (
                phase.UfEqn()
             ==
               *momentumTransfer[phase.name()]
              + fvOptions(alpha, rho, U)
            )
        );

        UEqns[phasei].relax();
        fvOptions.constrain(UEqns[phasei]);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
}
