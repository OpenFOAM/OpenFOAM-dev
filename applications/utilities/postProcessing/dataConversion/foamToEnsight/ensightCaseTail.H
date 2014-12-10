if (Pstream::master())
{
    ensightCaseFile << nl << "TIME" << nl
        << "time set:                      " << 1 << nl
        << "number of steps:               " << nTimeSteps << nl
        << "filename start number:         " << 0 << nl
        << "filename increment:            " << 1 << nl;

    ensightCaseFile << "time values:" << nl;

    ensightCaseFile.setf(ios_base::scientific, ios_base::floatfield);
    ensightCaseFile.precision(5);

    label count = 0;
    scalar Tcorr = 0.0;
    if (Times[0].value() < 0)
    {
        Tcorr = - Times[0].value();
        Info<< "Correcting time values. Adding " << Tcorr << endl;
    }

    forAll(Times, n)
    {
        ensightCaseFile << setw(12) << Times[n].value() + Tcorr << " ";

        if (++count % 6 == 0)
        {
            ensightCaseFile << nl;
        }
    }

    ensightCaseFile << nl;
}
