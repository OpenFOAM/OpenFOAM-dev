/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    stitchMesh

Description
    Utility to stitch or conform pairs of patches,
    converting the patch faces either into internal faces
    or conformal faces or another patch.

Usage
    \b stitchMesh (\<list of patch pairs\>)

    E.g. to stitch patches \c top1 to \c top2 and \c bottom1 to \c bottom2
        stitchMesh "((top1 top2) (bottom1 bottom2))"

    Options:
      - \par -overwrite \n
        Replace the old mesh with the new one, rather than writing the new one
        into a separate time directory

      - \par -region \<name\>
        Specify an alternative mesh region.

      - \par -fields
        Update vol and point fields

      - \par -tol
        Merge tolerance relative to local edge length (default 1e-4)

See also
    Foam::mergePatchPairs

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "timeSelector.H"
#include "ReadFields.H"
#include "volFields.H"
#include "pointFields.H"
#include "mergePatchPairs.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions(true, false);

    argList::addNote
    (
        "Stitch the faces on the specified patch pairs\n"
        "converting the overlapping patch faces into internal faces.\n"
    );

    argList::noParallel();
    #include "addOverwriteOption.H"
    #include "addMeshOption.H"
    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "fields",
        "update fields"
    );
    argList::addOption
    (
        "tol",
        "scalar",
        "merge tolerance relative to local edge length (default 1e-4)"
    );

    argList::validArgs.append("patchPairs");

    #include "setRootCase.H"
    #include "createTimeNoFunctionObjects.H"

    // Select time if specified
    timeSelector::selectIfPresent(runTime, args);

    #include "createSpecifiedMeshNoChangers.H"

    const scalar snapTol = args.optionLookupOrDefault("tol", 1e-4);
    const bool overwrite = args.optionFound("overwrite");
    const bool fields = args.optionFound("fields");

    const List<Pair<word>> patchPairNames((IStringStream(args[1])()));

    const word oldInstance = mesh.pointsInstance();

    // Search for list of objects for this time
    IOobjectList objects(mesh, runTime.name());

    if (fields) Info<< "Reading geometric fields" << nl << endl;

    #include "readVolFields.H"
    #include "readPointFields.H"

    Info<< endl;

    if (!overwrite)
    {
        runTime++;
    }

    // Stitch all the patch-pairs
    mergePatchPairs(mesh, patchPairNames, snapTol);

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    Info<< "Writing mesh to " << mesh.facesInstance() << endl;
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
