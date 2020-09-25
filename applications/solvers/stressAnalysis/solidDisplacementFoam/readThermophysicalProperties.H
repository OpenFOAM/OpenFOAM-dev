Info<< "Reading thermophysical properties\n" << endl;

solidDisplacementThermo thermo(mesh);

const volScalarField& E(thermo.E());
const volScalarField& nu(thermo.nu());

const volScalarField Cp(thermo.Cp());

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
