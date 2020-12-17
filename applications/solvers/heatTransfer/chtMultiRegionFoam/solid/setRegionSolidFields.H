const fvMesh& mesh = solidRegions[i];

solidThermo& thermo = thermoSolid[i];

tmp<volScalarField> trho = thermo.rho();
const volScalarField& rho = trho();

volScalarField& e = thermo.he();

const fv::options& fvOptions = fvOptionsSolid[i];

#include "checkRadiationModel.H"

solidNoLoopControl& pimple = pimples.solid(i);
