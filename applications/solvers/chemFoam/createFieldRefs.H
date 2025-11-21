basicChemistryModel& chemistry = pChemistry();
scalar dtChem = min(chemistry.deltaTChem()[0], runTime.deltaT().value());
PtrList<volScalarField>& Y = thermo.Y();
volScalarField& p = thermo.p();

const label nSpecie = thermo.Y().size();
const scalarList W(::W(thermo));

scalarList Y0(nSpecie, 0.0);

forAll(Y0, i)
{
    Y0[i] = thermo.Y(i)[0];
}

const scalar h0 = ::h0(thermo, Y0, p0, T0);
scalar rho0 = rho[0];
scalar u0 = h0 - p0/rho0;
scalar R0 = p0/(rho0*T0);
Rspecific[0] = R0;

scalar integratedHeat = 0.0;

Info << constProp << " will be held constant." << nl
    << " p   = " << p0 << " [Pa]" << nl
    << " T   = " << T0 << " [K] " << nl
    << " rho = " << rho0 << " [kg/m^3]" << nl
    << endl;
