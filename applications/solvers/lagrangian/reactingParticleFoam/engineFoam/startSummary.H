Info<< "Total cylinder mass: " << fvc::domainIntegrate(rho).value() << endl;
//Info<< "Total fuel mass: "
//    << fvc::domainIntegrate(rho*composition.ft()).value() << endl;

OFstream logSummaryFile
(
    runTime.path()/("logSummary." + runTime.timeName() + ".dat")
);

logSummaryFile
    << "# CA" << "      p" << "        T" << "       u'" << endl;
