Switch adjustTimeStep(runTime.controlDict().lookup("adjustTimeStep"));
scalar maxDeltaT(readScalar(runTime.controlDict().lookup("maxDeltaT")));
