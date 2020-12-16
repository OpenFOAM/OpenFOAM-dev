solidThermo& thermo = thermoSolid[i];

tmp<volScalarField> trho = thermo.rho();
const volScalarField& rho = trho();

volScalarField& e = thermo.he();

const fv::options& fvOptions = fvOptionsSolid[i];

solidNoLoopControl& pimple = pimples.solid(i);
