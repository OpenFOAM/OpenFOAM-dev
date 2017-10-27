{
    const scalar meanp = p.weightedAverage(mesh.V()).value();
    const scalar meanT = T.weightedAverage(mesh.V()).value();
    const scalar meanup =
        (sqrt((2.0/3.0)*turbulence->k()))().weightedAverage(mesh.V()).value();

    if (Pstream::master())
    {
        Info<< "Mean pressure:" << meanp << endl;
        Info<< "Mean temperature:" << meanT << endl;
        Info<< "Mean u':" << meanup << endl;

        logSummaryFile()
            << runTime.theta() << tab
            << meanp << tab
            << meanT << tab
            << meanup
            << endl;
    }
}
