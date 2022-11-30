{
    Info<< "k(" << runTime.name() << ") = "
        << 3.0/2.0*average(magSqr(U)).value() << endl;

    Info<< "epsilon(" << runTime.name() << ") = "
        << (0.5*nu*average(magSqr(twoSymm(fvc::grad(U))))).value() << endl;

    Info<< "U.f(" << runTime.name() << ") = "
        << 181.0*average(U & force).value() << endl;
}
