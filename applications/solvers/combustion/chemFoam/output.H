    runTime.write();

    Info<< "Sh = " << Sh
        << ", T = " << thermo.T()[0]
        << ", p = " << thermo.p()[0]
        << ", " << Y[0].name() << " = " << Y[0][0]
        << endl;

    post<< runTime.value() << token::TAB << thermo.T()[0] << token::TAB
        << thermo.p()[0] << endl;

