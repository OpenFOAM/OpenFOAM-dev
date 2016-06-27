IOdictionary propsDict
(
    IOobject
    (
        "particleTrackProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED
    )
);

const word cloudName(propsDict.lookup("cloudName"));

label sampleFrequency(readLabel(propsDict.lookup("sampleFrequency")));

label maxPositions(readLabel(propsDict.lookup("maxPositions")));

word setFormat(propsDict.lookupOrDefault<word>("setFormat", "vtk"));
