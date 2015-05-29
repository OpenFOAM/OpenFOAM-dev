    IOporosityModelList pZones(mesh);
    Switch pressureImplicitPorosity(false);

    // nUCorrectors used for pressureImplicitPorosity
    int nUCorr = 0;
    if (pZones.active())
    {
        // nUCorrectors for pressureImplicitPorosity
        nUCorr = simple.dict().lookupOrDefault<int>("nUCorrectors", 0);

        if (nUCorr > 0)
        {
            pressureImplicitPorosity = true;
            Info<< "Using pressure implicit porosity" << endl;
        }
        else
        {
            Info<< "Using pressure explicit porosity" << endl;
        }
    }
