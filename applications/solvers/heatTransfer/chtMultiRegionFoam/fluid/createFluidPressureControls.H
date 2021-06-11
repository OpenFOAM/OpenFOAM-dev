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

    hydrostaticInitialisation
    (
        p_rghFluid[i],
        thermoFluid[i].p(),
        rhoFluid[i],
        UFluid[i],
        ghFluid[i],
        ghfFluid[i],
        pRefFluid[i],
        thermoFluid[i],
        pimples.pimple(i).dict()
    );
}
