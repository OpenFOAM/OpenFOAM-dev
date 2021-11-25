    scalarField UMeanXvalues
    (
        channelIndexing.collapse(UMean.component(vector::X)())
    );

    scalarField UMeanYvalues
    (
        channelIndexing.collapse(UMean.component(vector::Y)())
    );

    scalarField UMeanZvalues
    (
        channelIndexing.collapse(UMean.component(vector::Z)())
    );

    scalarField RxxValues(channelIndexing.collapse(Rxx));
    scalarField RyyValues(channelIndexing.collapse(Ryy));
    scalarField RzzValues(channelIndexing.collapse(Rzz));
    scalarField RxyValues(channelIndexing.collapse(Rxy, true));

    scalarField pPrime2MeanValues(channelIndexing.collapse(pPrime2Mean));

    /*
    scalarField epsilonValues(channelIndexing.collapse(epsilonMean));

    scalarField nuMeanValues(channelIndexing.collapse(nuMean));
    scalarField nuPrimeValues(channelIndexing.collapse(nuPrime));

    scalarField gammaDotMeanValues(channelIndexing.collapse(gammaDotMean));
    scalarField gammaDotPrimeValues(channelIndexing.collapse(gammaDotPrime));
    */

    scalarField urmsValues(sqrt(mag(RxxValues)));
    scalarField vrmsValues(sqrt(mag(RyyValues)));
    scalarField wrmsValues(sqrt(mag(RzzValues)));

    scalarField kValues
    (
        0.5*(sqr(urmsValues) + sqr(vrmsValues) + sqr(wrmsValues))
    );

    setWriter::New(runTime.graphFormat())->write
    (
        runTime.globalPath()
       /functionObjects::writeFile::outputPrefix
       /args.executable()
       /runTime.timeName(),

       args.executable(),

       coordSet(true, "y", channelIndexing.y()),

       "Uf", UMeanXvalues,

       "u", urmsValues,
       "v", vrmsValues,
       "w", wrmsValues,
       "uv", RxyValues,

       "k", kValues,

       "pPrime2Mean", pPrime2MeanValues

       // "epsilon", epsilonValues,
       // "nu", nuMeanValues,
       // "nuPrime", nuPrimeValues,
       // "gammaDot", gammaDotMeanValues,
       // "gammaDotPrime", gammaDotPrimeValues
    );
