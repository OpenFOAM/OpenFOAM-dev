volScalarField::Internal contErr
(
    (fvc::ddt(rho) + fvc::div(mixture.rhoPhi()))()
);

forAll(mixture.phases(), phasei)
{
    const volScalarField& rho = mixture.phases()[phasei].thermo().rho();
    contErr -= (fvModels.source(mixture.phases()[phasei], rho)&rho);
}
