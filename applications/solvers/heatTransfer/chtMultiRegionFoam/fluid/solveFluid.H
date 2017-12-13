if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

if (frozenFlow)
{
    #include "EEqn.H"
}
else
{
    if (oCorr == 0)
    {
        #include "rhoEqn.H"
    }

    #include "UEqn.H"
    #include "YEqn.H"
    #include "EEqn.H"

    // --- PISO loop
    for (int corr=0; corr<nCorr; corr++)
    {
        #include "pEqn.H"
    }

    turbulence.correct();

    rho = thermo.rho();
}

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
