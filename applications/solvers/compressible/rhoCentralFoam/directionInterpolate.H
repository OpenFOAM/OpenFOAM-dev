namespace Foam
{

//- Interpolate field vf according to direction dir
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& dir,
    const word& reconFieldName = word::null
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf
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

    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();

    sf.rename(vf.name() + '_' + dir.name());

    return tsf;
}

}
