const volScalarField& psi = thermo.psi();
const volScalarField& T = thermo.T();
filmModelType& surfaceFilm = tsurfaceFilm();
const label inertIndex(composition.species()[inertSpecie]);
