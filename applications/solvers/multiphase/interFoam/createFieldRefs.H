volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

twoPhaseChangeModel& phaseChange = phaseChangePtr();
tmp<volScalarField> rAU;
