/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    surfaceSplitByPatch

Description
    Writes regions of triSurface to separate files.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "triSurface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "write surface mesh regions to separate files"
    );

    argList::noParallel();
    argList::validArgs.append("input surfaceFile");
    argList args(argc, argv);

    const fileName surfName = args[1];

    Info<< "Reading surf from " << surfName << " ..." << nl << endl;

    fileName surfBase = surfName.lessExt();

    word extension = surfName.ext();

    triSurface surf(surfName);

    Info<< "Writing regions to separate files ..."
        << nl << endl;


    const geometricSurfacePatchList& patches = surf.patches();

    forAll(patches, patchI)
    {
        const geometricSurfacePatch& pp = patches[patchI];

        word patchName = pp.name();

        if (patchName.empty())
        {
            patchName = "patch" + Foam::name(patchI);
        }

        fileName outFile(surfBase + '_' + patchName + '.' + extension);

        Info<< "   Writing patch " << patchName << " to file " << outFile
            << endl;


        // Collect faces of region
        boolList includeMap(surf.size(), false);

        forAll(surf, faceI)
        {
            const labelledTri& f = surf[faceI];

            if (f.region() == patchI)
            {
                includeMap[faceI] = true;
            }
        }

        // Subset triSurface
        labelList pointMap;
        labelList faceMap;

        triSurface subSurf
        (
            surf.subsetMesh
            (
                includeMap,
                pointMap,
                faceMap
            )
        );

        subSurf.write(outFile);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
