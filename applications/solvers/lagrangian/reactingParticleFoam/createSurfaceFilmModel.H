Info<< "\nConstructing surface film model" << endl;

autoPtr<regionModels::surfaceFilmModel> tsurfaceFilm
(
    regionModels::surfaceFilmModel::New(mesh, g)
);
