namespace Foam
{

//- Interpolate field vf according to direction dir
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& dir,
    const word& reconFieldName = word::null
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        fvc::interpolate
        (
            vf,
            dir,
            "reconstruct("
          + (reconFieldName != word::null ? reconFieldName : vf.name())
          + ')'
        )
    );

    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf();

    sf.rename(vf.name() + '_' + dir.name());

    // Correct BCs of the positive (outgoing) fluxes
    forAll(sf.boundaryField(), patchi)
    {
        if
        (
            !sf.boundaryField()[patchi].coupled()
         && sf.boundaryField()[patchi].size()
         && !vf.boundaryField()[patchi].fixesValue()
         && dir.boundaryField()[patchi][0] > 0
        )
        {
            sf.boundaryField()[patchi] =
                vf.boundaryField()[patchi].patchInternalField();
        }
    }

    return tsf;
}


template<class Type>
wordList derivedPatchFieldTypes
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    wordList phiTypes
    (
        vf.boundaryField().size(),
        calculatedFvPatchField<Type>::typeName
    );

    forAll(vf.boundaryField(), patchi)
    {
        if (vf.boundaryField()[patchi].fixesValue())
        {
            phiTypes[patchi] = fixedValueFvPatchField<Type>::typeName;
        }
    }

    return phiTypes;
}

}
