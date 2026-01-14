const Foam::wordList regionNames
(
    args.optionFound("allRegions")
  ? runTime.regionNames()
  : wordList
    (
        1,
        args.optionFound("region")
      ? args.optionRead<word>("region")
      : polyMesh::defaultRegion
    )
);
