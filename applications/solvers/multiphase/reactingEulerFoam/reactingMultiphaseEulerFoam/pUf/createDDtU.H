    PtrList<surfaceScalarField> ddtPhis(phases.size());

    forAll(phases, phasei)
    {

        phaseModel& phase = phases[phasei];

        ddtPhis.set
        (
            phasei,
            (
                (phase.phi() - phase.phi().oldTime())/runTime.deltaT()
            ).ptr()
        );
    }
