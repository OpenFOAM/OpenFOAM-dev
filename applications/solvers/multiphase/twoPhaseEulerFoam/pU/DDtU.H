{
    DDtU1 =
        fvc::ddt(U1)
      + fvc::div(phi1, U1)
      - fvc::div(phi1)*U1;
    mrfZones.addCoriolis(U1, DDtU1);

    DDtU2 =
        fvc::ddt(U2)
      + fvc::div(phi2, U2)
      - fvc::div(phi2)*U2;
    mrfZones.addCoriolis(U2, DDtU2);
}
