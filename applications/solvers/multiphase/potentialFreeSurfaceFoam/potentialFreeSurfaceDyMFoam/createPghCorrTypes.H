    wordList pcorrTypes
    (
        p_gh.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    for (label i=0; i<p_gh.boundaryField().size(); i++)
    {
        if (p_gh.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }
