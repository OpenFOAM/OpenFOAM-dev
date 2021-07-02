IOdictionary propsDict(systemDict("particleTracksDict", args, runTime));

const word cloudName(propsDict.lookup("cloudName"));

label sampleFrequency(propsDict.lookup<label>("sampleFrequency"));

label maxPositions(propsDict.lookup<label>("maxPositions"));

word setFormat(propsDict.lookupOrDefault<word>("setFormat", "vtk"));
