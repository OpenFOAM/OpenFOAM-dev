    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];

        ddtPhis[phasei] =
            (phase.phi() - phase.phi().oldTime())/runTime.deltaT();
    }
