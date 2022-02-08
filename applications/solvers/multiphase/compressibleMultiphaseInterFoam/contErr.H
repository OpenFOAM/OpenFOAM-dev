volScalarField::Internal contErr
(
    (fvc::ddt(rho) + fvc::div(mixture.rhoPhi()))()
);

forAllConstIter
(
    PtrDictionary<phaseModel>,
    mixture.phases(),
    phase
)
{
    const volScalarField& rho = phase().thermo().rho()();
    contErr -= (fvModels.source(phase(), rho)&rho);
}
