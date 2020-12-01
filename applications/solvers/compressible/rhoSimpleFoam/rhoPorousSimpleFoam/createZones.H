    IOporosityModelList pZones(mesh);
    Switch pressureImplicitPorosity(false);

    // nUCorrectors used for pressureImplicitPorosity
    int nUCorr = 0;
    if (pZones.UPtrList<porosityModel>::size())
    {
        // nUCorrectors for pressureImplicitPorosity
        simple.dict().readIfPresent("nUCorrectors", nUCorr);

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
