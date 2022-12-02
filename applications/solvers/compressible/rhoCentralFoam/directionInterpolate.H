namespace Foam
{

//- Interpolate field vf according to direction dir
template<class Type>
tmp<SurfaceField<Type>> interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& dir,
    const word& reconFieldName = word::null
)
{
    tmp<SurfaceField<Type>> tsf
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

    SurfaceField<Type>& sf = tsf.ref();

    sf.rename(vf.name() + '_' + dir.name());

    return tsf;
}

}
