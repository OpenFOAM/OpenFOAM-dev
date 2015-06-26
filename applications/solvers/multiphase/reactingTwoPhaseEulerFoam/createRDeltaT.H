bool LTS =
    word(mesh.ddtScheme("default"))
 == fv::localEulerDdtScheme<scalar>::typeName;

tmp<volScalarField> trDeltaT;
tmp<volScalarField> trSubDeltaT;

if (LTS)
{
    scalar maxDeltaT
    (
        pimple.dict().lookupOrDefault<scalar>("maxDeltaT", GREAT)
    );

    trDeltaT = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rDeltaT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    trSubDeltaT = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rSubDeltaT",
                runTime.timeName(),
                mesh
            ),
            mesh,
            1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT)
        )
    );
}
