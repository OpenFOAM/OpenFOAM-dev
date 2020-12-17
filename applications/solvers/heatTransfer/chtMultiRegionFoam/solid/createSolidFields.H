// Initialise solid field pointer lists
PtrList<solidThermo> thermoSolid(solidRegions.size());
PtrList<fv::options> fvOptionsSolid(solidRegions.size());

// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "*** Reading solid mesh thermophysical properties for region "
        << solidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermoSolid\n" << endl;
    thermoSolid.set(i, solidThermo::New(solidRegions[i]));

    Info<< "    Adding to fvOptionsSolid\n" << endl;
    fvOptionsSolid.set
    (
        i,
        new fv::options(solidRegions[i])
    );
}
