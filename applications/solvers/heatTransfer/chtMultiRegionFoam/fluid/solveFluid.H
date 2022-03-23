if (!pimple.flow())
{
    if (pimple.models())
    {
        fvModels.correct();
    }

    if (pimple.thermophysics())
    {
        tmp<fv::convectionScheme<scalar>> mvConvection(nullptr);

        if (Ecorr == 0)
        {
            #include "YEqn.H"
        }
        #include "EEqn.H"
    }
}
else
{
    if (Ecorr == 0)
    {
        if (!mesh.schemes().steady() && pimples.firstPimpleIter())
        {
            #include "rhoEqn.H"
        }

        if (pimple.models())
        {
            fvModels.correct();
        }

        #include "UEqn.H"
    }

    if (pimple.thermophysics())
    {
        tmp<fv::convectionScheme<scalar>> mvConvection(nullptr);

        if (Ecorr == 0)
        {
            #include "YEqn.H"
        }
        #include "EEqn.H"
    }

    if (Ecorr == nEcorr - 1)
    {
        tmp<fvVectorMatrix>& tUEqn = UEqns[i];
        fvVectorMatrix& UEqn = tUEqn.ref();

        // --- PISO loop
        while (pimple.correct())
        {
            #include "../../buoyantFoam/pEqn.H"
        }

        if (pimples.pimpleTurbCorr(i))
        {
            turbulence.correct();
            thermophysicalTransport.correct();
        }

        if (!mesh.schemes().steady() && pimples.finalPimpleIter())
        {
            rho = thermo.rho();
        }
    }
}
