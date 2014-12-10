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

    const label nIntervals(readLabel(pdfDictionary.lookup("nIntervals")));

    const label nSamples(readLabel(pdfDictionary.lookup("nSamples")));

    const bool writeData(readBool(pdfDictionary.lookup("writeData")));


    const fileName pdfPath = runTime.path()/"pdf";
    mkDir(pdfPath);

    cachedRandom rndGen(label(0), -1);

    autoPtr<distributionModels::distributionModel> p
    (
        distributionModels::distributionModel::New
        (
            pdfDictionary,
            rndGen
        )
    );

    const scalar xMin = p->minValue();
    const scalar xMax = p->maxValue();

    autoPtr<OFstream> filePtr(NULL);
    if (writeData)
    {
        fileName fName = pdfPath/(p->type() + ".data");
        Info<< "Writing " << p->type() << " data samples to file:" << nl
            << fName << nl << endl;

        filePtr.reset(new OFstream(fName));
    }

    scalarField samples(nIntervals, 0);

