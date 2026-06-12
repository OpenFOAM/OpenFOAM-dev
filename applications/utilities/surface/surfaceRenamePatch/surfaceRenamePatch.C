/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    surfaceRenamePatch

Description
    Rename a patch in a surface geometry file. The patch can be identified by
    name or a point close to the patch.

Usage
    \b surfaceRenamePatch name inputFile outputFile [OPTION]

    Options:
      - \par -point
        Point close to the patch

      - \par -name
        Current name of the patch

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "MeshedSurfaces.H"
#include "indexedOctree.H"
#include "treeDataPrimitivePatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    argList::validArgs.append("new name");
    argList::validArgs.append("input surface file");
    argList::validArgs.append("output surface file");

    argList::addOption
    (
        "point",
        "point",
        "point close to the patch; e.g., '(0 0 0.5)'"
    );

    argList::addOption
    (
        "name",
        "word",
        "current name of the patch; e.g., patch0"
    );

    argList::addNote
    (
        "Rename an existing patch to <new name> "
        "in a surface geometry file.\n\n"
        "The existing patch can be identified by one of two options:\n"
        "  -name <word>: selects the patch by name\n"
        "  -point <point>: select the nearest patch to the specified point\n\n"
        "Examples:\n"
        "+ To rename 'patch0' to 'inlet'\n"
        "  surfaceRenamePatch -name patch0 inlet in.obj out.obj\n"
        "+ To rename the patch nearest to the point '(1 2 3)' to 'inlet'\n"
        "  surfaceRenamePatch -point '(1 2 3)' inlet in.obj out.obj\n"
    );

    argList args(argc, argv);

    const word newPatchName(args.argRead<word>(1));
    const fileName inFile(args.argRead<fileName>(2));
    const fileName outFile(args.argRead<fileName>(3));

    if (!args.optionFound("point") && !args.optionFound("name"))
    {
        FatalErrorInFunction
            << "Missing option: -name or -point"
            << exit(FatalError);
    }

    if (args.optionFound("point") && args.optionFound("name"))
    {
        FatalErrorInFunction
            << "Both options provided: -patch and -point"
            << exit(FatalError);
    }

    Info<< "Reading surface " << inFile << nl << endl;

    meshedSurface surf(inFile);

    // Get the index of the zone ...
    label surfZonei = -1;
    if (args.optionFound("name"))
    {
        const word surfZoneName = args.optionRead<word>("name");

        // Loop through all the zones looking for the one with the given name
        forAll(surf.surfZones(), i)
        {
            const surfZone& sZone = surf.surfZones()[i];

            if (sZone.name() == surfZoneName)
            {
                surfZonei = i;

                break;
            }
        }

        // Error if the zone was not found
        if (surfZonei == -1)
        {
            FatalErrorInFunction
                << "Patch " << surfZoneName << " not found"
                << exit(FatalError);
        }
    }
    else // (args.optionFound("point"))
    {
        const point p = args.optionRead<point>("point");

        // Construct a search tree for the surface
        typedef treeDataPrimitivePatch<meshedSurface> treeType;
        const indexedOctree<treeType> tree
        (
            treeType
            (
                false,
                surf,
                indexedOctree<treeType>::perturbTol()
            ),
            treeBoundBox(surf.points()).extend(1e-4),
            8,
            10,
            3
        );

        // Find the nearest face to the tree
        const pointIndexHit hit = tree.findNearest(p, sqr(great));

        if (!hit.hit())
        {
            FatalErrorInFunction
                << "Patch nearest to " << p << " not found"
                << exit(FatalError);
        }

        // Loop through all the zones looking for the one which contains the
        // nearest face
        forAll(surf.surfZones(), i)
        {
            const surfZone& sZone = surf.surfZones()[i];

            if
            (
                hit.index() >= sZone.start()
             && hit.index() < sZone.start() + sZone.size()
            )
            {
                surfZonei = i;

                Info<< "Found patch " << sZone.name()
                    << " nearest to " << p << nl << endl;

                break;
            }
        }
    }

    // Rename the patch
    Info<< "Renaming '" << surf.surfZones()[surfZonei].name()
        << "' to '" << newPatchName << "'"
        << nl << endl;

    surf.surfZones()[surfZonei].name() = newPatchName;

    Info<< "Writing surface " << outFile << nl << endl;
    surf.write(outFile);

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
