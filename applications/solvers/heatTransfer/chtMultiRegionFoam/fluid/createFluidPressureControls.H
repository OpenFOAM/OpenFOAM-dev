PtrList<pressureReference> pressureReferenceFluid(fluidRegions.size());

forAll(fluidRegions, i)
{
    pressureReferenceFluid.set
    (
        i,
        new pressureReference
        (
            p_rghFluid[i],
            pimples.pimple(i).dict(),
            false
        )
    );
}
