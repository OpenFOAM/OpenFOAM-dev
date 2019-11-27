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

label sampleFrequency(propsDict.lookup<label>("sampleFrequency"));

label maxPositions(propsDict.lookup<label>("maxPositions"));

word setFormat(propsDict.lookupOrDefault<word>("setFormat", "vtk"));
