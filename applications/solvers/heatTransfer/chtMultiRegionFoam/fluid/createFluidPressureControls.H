PtrList<pressureControl> pressureControlFluid(fluidRegions.size());

forAll(fluidRegions, i)
{
    pressureControlFluid.set
    (
        i,
        new pressureControl
        (
            p_rghFluid[i],
            rhoFluid[i],
            pimples.pimple(i).dict(),
            false
        )
    );
}
