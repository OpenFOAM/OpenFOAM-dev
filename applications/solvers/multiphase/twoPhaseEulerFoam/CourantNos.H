#   include "CourantNo.H"

{
    scalar UrCoNum = 0.5*gMax
    (
        fvc::surfaceSum(mag(phi1 - phi2))().internalField()/mesh.V().field()
    )*runTime.deltaTValue();

    Info<< "Max Ur Courant Number = " << UrCoNum << endl;

    CoNum = max(CoNum, UrCoNum);
}
