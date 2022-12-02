namespace Foam
{

//- Interpolate field vf according to direction dir
template<class Type>
tmp<SurfaceField<Type>> interpolate
(
    const VolField<Type>& vf,
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
