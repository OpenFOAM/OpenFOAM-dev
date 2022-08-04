const Foam::wordList regionNames
(
    args.optionFound("allRegions")
  ? runTime .controlDict().subDict("regionSolvers").toc()
  : wordList
    (
        1,
        args.optionFound("region")
      ? args.optionRead<word>("region")
      : polyMesh::defaultRegion
    )
);
