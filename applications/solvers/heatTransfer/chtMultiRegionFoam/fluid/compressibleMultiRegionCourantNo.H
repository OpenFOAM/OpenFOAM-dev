    scalar CoNum = -GREAT;

    forAll(fluidRegions, regionI)
    {
        CoNum = max
        (
            compressibleCourantNo
            (
                fluidRegions[regionI],
                runTime,
                rhoFluid[regionI],
                phiFluid[regionI]
            ),
            CoNum
        );
    }
/*
    forAll (porousFluidRegions, porousI)
    {
        CoNum = max
        (
            compressibleCourantNo
            (
                porousFluidRegions[porousI],
                runTime,
                rhoPorous[porousI],
                phiPorous[porousI]
            ),
            CoNum
        );
    }
*/