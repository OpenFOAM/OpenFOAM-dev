label inertIndex = -1;
if (Y.size())
{
    const word inertSpecie(thermo.lookup("inertSpecie"));
    if (!composition.species().found(inertSpecie))
    {
        FatalIOErrorIn(args.executable().c_str(), thermo)
            << "Inert specie " << inertSpecie
            << " not found in available species "
            << composition.species()
            << exit(FatalIOError);
    }
    inertIndex = composition.species()[inertSpecie];
}

const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();
regionModels::surfaceFilmModel& surfaceFilm = tsurfaceFilm();

Switch solvePrimaryRegion
(
    pimple.dict().lookupOrDefault<Switch>("solvePrimaryRegion", true)
);
