if (!pimple.solveFlow())
{
    #include "EEqn.H"
}
else
{
    if (!mesh.steady() && pimples.nCorrPIMPLE() <= 1)
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

    if (!mesh.steady() && pimples.finalIter())
    {
        rho = thermo.rho();
    }
}
