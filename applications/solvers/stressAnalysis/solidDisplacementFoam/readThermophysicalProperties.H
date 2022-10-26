Info<< "Reading thermophysical properties\n" << endl;

solidDisplacementThermo thermo(mesh);

autoPtr<solidThermophysicalTransportModel> thermophysicalTransport;

if (thermo.thermalStress())
{
    thermophysicalTransport = solidThermophysicalTransportModel::New(thermo);
}

const volScalarField& E(thermo.E());
const volScalarField& nu(thermo.nu());

Info<< "Calculating Lame's coefficients\n" << endl;

const volScalarField mu(E/(2*(1 + nu)));

const volScalarField lambda
(
    thermo.planeStress()
  ? nu*E/((1 + nu)*(1 - nu))
  : nu*E/((1 + nu)*(1 - 2*nu))
);

const volScalarField threeK
(
    thermo.planeStress()
  ? E/(1 - nu)
  : E/(1 - 2*nu)
);

const volScalarField threeKalpha("threeKalpha", threeK*thermo.alphav());
