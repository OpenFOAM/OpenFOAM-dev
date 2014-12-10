Info<< "\nConstructing surface film model" << endl;

typedef regionModels::surfaceFilmModels::surfaceFilmModel filmModelType;

autoPtr<filmModelType> tsurfaceFilm(filmModelType::New(mesh, g));
filmModelType& surfaceFilm = tsurfaceFilm();

