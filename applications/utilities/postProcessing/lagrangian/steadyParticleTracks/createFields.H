const word dictName("particleTrackDict");

#include "setConstantMeshDictionaryIO.H"

IOdictionary propsDict(dictIO);

word cloudName(propsDict.lookup("cloudName"));

List<word> userFields(propsDict.lookup("fields"));
