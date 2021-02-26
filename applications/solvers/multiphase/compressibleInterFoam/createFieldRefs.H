twoPhaseChangeModel& phaseChange = phaseChangePtr();

volScalarField& alpha2(mixture.alpha2());

const volScalarField& rho1 = mixture.thermo1().rho();
const volScalarField& rho2 = mixture.thermo2().rho();

volScalarField& p = mixture.p();
volScalarField& T = mixture.T();
const volScalarField& psi1 = mixture.thermo1().psi();
const volScalarField& psi2 = mixture.thermo2().psi();

tmp<volScalarField> rAU;

if (correctPhi)
{
    rAU = new volScalarField
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimTime/dimDensity, 1)
    );
}
