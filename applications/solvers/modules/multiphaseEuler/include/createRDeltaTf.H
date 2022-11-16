tmp<surfaceScalarField> trDeltaTf;

if (LTS && faceMomentum)
{
    trDeltaTf = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                fv::localEulerDdt::rDeltaTfName,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, 1)
        )
    );
}
