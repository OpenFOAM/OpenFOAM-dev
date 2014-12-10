Info<< "Mean pressure:" << p.weightedAverage(mesh.V()).value() << endl;
Info<< "Mean temperature:" << thermo.T().weightedAverage(mesh.V()).value()
    << endl;
Info<< "Mean u':"
    << (sqrt((2.0/3.0)*turbulence->k()))().weightedAverage(mesh.V()).value()
    << endl;

logSummaryFile
    << runTime.theta() << tab
    << p.weightedAverage(mesh.V()).value() << tab
    << thermo.T().weightedAverage(mesh.V()).value() << tab
    << (sqrt((2.0/3.0)*turbulence->k()))().weightedAverage(mesh.V()).value()
    << tab
    << 1 - b.weightedAverage(mesh.V()).value()
    << endl;
