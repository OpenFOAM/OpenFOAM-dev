Info<< "\nConstructing surface film model" << endl;

typedef regionModels::surfaceFilmModels::surfaceFilmModel filmModelType;

SLGThermo slgThermo(mesh, mixture.thermo1());
autoPtr<filmModelType> tsurfaceFilm(filmModelType::New(mesh, g));
