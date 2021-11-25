    IOdictionary pdfDictionary
    (
        IOobject
        (
            "pdfDict",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const label nIntervals(pdfDictionary.lookup<label>("nIntervals"));

    const label nSamples(pdfDictionary.lookup<label>("nSamples"));

    const bool writeData(readBool(pdfDictionary.lookup("writeData")));

    const fileName pdfPath =
        runTime.globalPath()
       /functionObjects::writeFile::outputPrefix
       /args.executable()
       /runTime.timeName();

    mkDir(pdfPath);

    Random rndGen(label(0));

    autoPtr<distributionModel> p
    (
        distributionModel::New
        (
            pdfDictionary,
            rndGen
        )
    );

    const scalar xMin = p->minValue();
    const scalar xMax = p->maxValue();

    autoPtr<OFstream> filePtr(nullptr);

    if (writeData)
    {
        fileName fName = pdfPath/(p->type() + ".data");
        Info<< "Writing " << p->type() << " data samples to file:" << nl
            << fName << nl << endl;

        filePtr.reset(new OFstream(fName));
    }

    scalarField samples(nIntervals, 0);
