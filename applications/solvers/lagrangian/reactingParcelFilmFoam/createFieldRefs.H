const label inertIndex(composition.species()[inertSpecie]);

const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();
filmModelType& surfaceFilm = tsurfaceFilm();
