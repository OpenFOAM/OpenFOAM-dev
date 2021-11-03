Info<< "Mean pressure:" << p.weightedAverage(mesh.V()).value() << endl;
Info<< "Mean temperature:" << T.weightedAverage(mesh.V()).value() << endl;
Info<< "Mean u':"
    << (sqrt((2.0/3.0)*turbulence->k()))().weightedAverage(mesh.V()).value()
    << endl;

logSummaryFile
    << runTime.userTimeValue() << tab
    << p.weightedAverage(mesh.V()).value() << tab
    << T.weightedAverage(mesh.V()).value() << tab
    << (sqrt((2.0/3.0)*turbulence->k()))().weightedAverage(mesh.V()).value()
    << endl;
