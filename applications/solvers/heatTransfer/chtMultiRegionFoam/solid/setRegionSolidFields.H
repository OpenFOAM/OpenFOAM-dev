    fvMesh& mesh = solidRegions[i];
    solidThermo& thermo = thermos[i];
    const radiation::radiationModel& radiation = radiations[i];

    tmp<volScalarField> trho = thermo.rho();
    const volScalarField& rho = trho();

    tmp<volScalarField> tcp = thermo.Cp();
    const volScalarField& cp = tcp();

    tmp<volSymmTensorField> taniAlpha;
    if (!thermo.isotropic())
    {
        volSymmTensorField& aniAlpha = aniAlphas[i];
        tmp<volVectorField> tkappaByCp = thermo.Kappa()/cp;
        const coordinateSystem& coodSys = coordinates[i];

        aniAlpha.internalField() =
            coodSys.R().transformVector(tkappaByCp());
        aniAlpha.correctBoundaryConditions();

        taniAlpha = tmp<volSymmTensorField>
        (
            new volSymmTensorField(aniAlpha)
        );
    }


    volScalarField& h = thermo.he();

    const volScalarField& betav = betavSolid[i];

    fv::IOoptionList& fvOptions = solidHeatSources[i];
