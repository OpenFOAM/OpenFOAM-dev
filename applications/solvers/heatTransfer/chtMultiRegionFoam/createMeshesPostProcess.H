#include "createMeshes.H"

if (!fluidRegions.size())
{
    FatalErrorIn(args.executable())
        << "No fluid meshes present" << exit(FatalError);
}

fvMesh& mesh = fluidRegions[0];
