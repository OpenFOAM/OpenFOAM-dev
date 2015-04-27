    Info<< "Calculating field DDtU1 and DDtU2\n" << endl;

    volVectorField DDtU1
    (
        "DDtU1",
        fvc::ddt(U1)
      + fvc::div(phi1, U1)
      - fvc::div(phi1)*U1
    );

    volVectorField DDtU2
    (
        "DDtU2",
        fvc::ddt(U2)
      + fvc::div(phi2, U2)
      - fvc::div(phi2)*U2
    );
