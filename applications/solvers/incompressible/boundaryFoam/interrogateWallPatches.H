// Search for wall patches faces and store normals

label faceId(-1);
label patchId(-1);
label nWallFaces(0);
vector wallNormal(vector::zero);

const fvPatchList& patches = mesh.boundary();

forAll(patches, patchi)
{
    const fvPatch& currPatch = patches[patchi];

    if (isA<wallFvPatch>(currPatch))
    {
        const vectorField nf(currPatch.nf());

        forAll(nf, facei)
        {
            nWallFaces++;

            if (nWallFaces == 1)
            {
                wallNormal = -nf[facei];
                faceId = facei;
                patchId = patchi;
            }
            else if (nWallFaces == 2)
            {
                const vector wallNormal2 = -nf[facei];

                //- Check that wall faces are parallel
                if
                (
                    mag(wallNormal & wallNormal2) > 1.01
                 || mag(wallNormal & wallNormal2) < 0.99
                )
                {
                    FatalErrorIn(args.executable())
                        << "wall faces are not parallel for patches "
                        << patches[patchId].name() << " and "
                        << currPatch.name() << nl
                        << exit(FatalError);
                }
            }
            else
            {
                FatalErrorIn(args.executable()) << "number of wall faces > 2"
                    << nl << exit(FatalError);
            }
        }
    }
}

if (nWallFaces == 0)
{
    FatalErrorIn(args.executable()) << "No wall patches identified"
        << exit(FatalError);
}
else
{
    Info<< "Generating wall data for patch: " << patches[patchId].name()
        << endl;
}

// store local id of near-walll cell to process
label cellId = patches[patchId].faceCells()[faceId];

// create position array for graph generation
scalarField y
(
    wallNormal
  & (mesh.C().internalField() - mesh.C().boundaryField()[patchId][faceId])
);

Info<< "    Height to first cell centre y0 = " << y[cellId] << endl;
