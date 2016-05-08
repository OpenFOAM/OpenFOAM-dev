volScalarField& alpha2 = phase2;

surfaceScalarField& phi1 = phase1.phi();
surfaceScalarField& alphaPhi1 = phase1.alphaPhi();
surfaceScalarField& alphaRhoPhi1 = phase1.alphaRhoPhi();

surfaceScalarField& phi2 = phase2.phi();
surfaceScalarField& alphaPhi2 = phase2.alphaPhi();
surfaceScalarField& alphaRhoPhi2 = phase2.alphaRhoPhi();

surfaceScalarField& phi = fluid.phi();

rhoThermo& thermo1 = phase1.thermo();
rhoThermo& thermo2 = phase2.thermo();

volScalarField& rho1 = thermo1.rho();
const volScalarField& psi1 = thermo1.psi();

volScalarField& rho2 = thermo2.rho();
const volScalarField& psi2 = thermo2.psi();
