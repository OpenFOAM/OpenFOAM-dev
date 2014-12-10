    autoPtr<radiation::radiationModel> radiation
    (
        radiation::radiationModel::New(thermo.T())
    );
