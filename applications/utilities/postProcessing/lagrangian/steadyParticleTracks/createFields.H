IOdictionary propsDict(systemDict("steadyParticleTracksDict", args, runTime));

word cloudName(propsDict.lookup("cloudName"));

List<word> userFields(propsDict.lookup("fields"));
