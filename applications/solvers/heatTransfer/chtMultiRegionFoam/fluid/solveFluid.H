if (pimple.frozenFlow())
{
    fvModels.correct();

    #include "YEqn.H"
    #include "EEqn.H"
}
else
{
    if (!mesh.steady() && pimples.firstPimpleIter())
    {
        #include "rhoEqn.H"
    }

    fvModels.correct();

    #include "UEqn.H"
    #include "YEqn.H"
    #include "EEqn.H"

    // --- PISO loop
    while (pimple.correct())
    {
        #include "../../buoyantPimpleFoam/pEqn.H"
    }

    if (pimples.pimpleTurbCorr(i))
    {
        turbulence.correct();
        thermophysicalTransport.correct();
    }

    if (!mesh.steady() && pimples.finalPimpleIter())
    {
        rho = thermo.rho();
    }
}
