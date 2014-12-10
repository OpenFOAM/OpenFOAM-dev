    wordList pcorrTypes
    (
        p_rgh.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    for (label i=0; i<p_rgh.boundaryField().size(); i++)
    {
        if (p_rgh.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }
