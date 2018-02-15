Info<< "Constructing face momentum equations" << endl;

PtrList<fvVectorMatrix> UEqns(phases.size());

{
    // Update of coefficients
    fluid.momentumTransfer();

    PtrList<fvVectorMatrix> UgradUs(phases.size());

    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];

        UgradUs.set
        (
            phasei,
            new fvVectorMatrix
            (
                fvm::div(phase.phi(), phase.U())
              - fvm::Sp(fvc::div(phase.phi()), phase.U())
              + MRF.DDt(phase.U())
            )
        );
    }

    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];

        const volScalarField& alpha = phase;
        volScalarField& rho = phase.thermo().rho();
        const surfaceScalarField& alphaRhoPhi = phase.alphaRhoPhi();
        volVectorField& U = phase.U();


        UEqns.set
        (
            phasei,
            new fvVectorMatrix
            (
                fvm::div(alphaRhoPhi, U) - fvm::Sp(fvc::div(alphaRhoPhi), U)
              + MRF.DDt(alpha*rho, U)
              + phase.divDevRhoReff()
              + fvm::SuSp(fvOptions(alpha, rho)&rho, U)
             ==
                fvOptions(alpha, rho, U)
            )
        );

        // Add gradient part of virtual mass force
        forAllConstIter
        (
            phaseSystem::VmTable,
            fluid.Vms(),
            VmIter
        )
        {
            const volScalarField& Vm(*VmIter());

            const phasePair& pair(fluid.phasePairs()[VmIter.key()]);

            if (pair.contains(phase))
            {
                const phaseModel& otherPhase = pair.otherPhase(phase);

                UEqns[phasei] +=
                    Vm
                   *(
                        UgradUs[phasei]
                      - (UgradUs[otherPhase.index()] & otherPhase.U())
                    );
            }
        }

        UEqns[phasei].relax();
        fvOptions.constrain(UEqns[phasei]);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
}
