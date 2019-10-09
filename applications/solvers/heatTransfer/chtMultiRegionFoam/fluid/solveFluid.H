if (pimple.frozenFlow())
{
    #include "EEqn.H"
}
else
{
    if (!mesh.steady() && pimples.firstPimpleIter())
    {
        #include "rhoEqn.H"
    }

    #include "UEqn.H"
    #include "YEqn.H"
    #include "EEqn.H"

    // --- PISO loop
    while (pimple.correct())
    {
        #include "pEqn.H"
    }

    if (pimples.pimpleTurbCorr(i))
    {
        turbulence.correct();
    }

    if (!mesh.steady() && pimples.finalPimpleIter())
    {
        rho = thermo.rho();
    }
}
