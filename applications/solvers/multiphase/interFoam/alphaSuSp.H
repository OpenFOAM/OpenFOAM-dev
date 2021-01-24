// Phase change alpha1 source
Pair<tmp<volScalarField::Internal>> phaseChangeS
(
    phaseChange.Salpha(alpha1)
);
const tmp<volScalarField::Internal>& Su = phaseChangeS[0];
const tmp<volScalarField::Internal>& Sp = phaseChangeS[1];

tmp<volScalarField> divU;

if (Sp.valid())
{
    divU =
    (
        mesh.moving()
      ? fvc::div(phiCN() + mesh.phi())
      : fvc::div(phiCN())
    );
}
