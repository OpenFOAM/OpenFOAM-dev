Info<< "Total cylinder mass: " << fvc::domainIntegrate(rho).value() << endl;

autoPtr<OFstream> logSummaryFile;

if (Pstream::master())
{
    logSummaryFile = new OFstream
    (
        runTime.rootPath()/runTime.globalCaseName()
       /("logSummary." + runTime.timeName() + ".dat")
    );

    logSummaryFile()
        << "# CA" << "      p" << "        T" << "       u'" << "        c"
        << endl;
}
