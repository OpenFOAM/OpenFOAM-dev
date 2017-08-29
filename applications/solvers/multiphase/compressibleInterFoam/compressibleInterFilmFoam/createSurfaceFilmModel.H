Info<< "\nConstructing surface film model" << endl;

SLGThermo slgThermo(mesh, mixture.thermo1());
autoPtr<regionModels::surfaceFilmModel> tsurfaceFilm
(
    regionModels::surfaceFilmModel::New(mesh, g)
);
