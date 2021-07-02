const word dictName("particleTrackDict");

IOdictionary propsDict
(
    systemDict(dictName, args, runTime, word::null, runTime.constant())
);

word cloudName(propsDict.lookup("cloudName"));

List<word> userFields(propsDict.lookup("fields"));
