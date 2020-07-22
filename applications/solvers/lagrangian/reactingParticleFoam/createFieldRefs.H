const label inertIndex(composition.species()[inertSpecie]);

const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();
regionModels::surfaceFilmModel& surfaceFilm = tsurfaceFilm();

Switch solvePrimaryRegion
(
    pimple.dict().lookupOrDefault<Switch>("solvePrimaryRegion", true)
);
