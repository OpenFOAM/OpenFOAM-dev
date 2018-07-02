/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    surfaceAdd

Description
    Add two surfaces. Does geometric merge on points. Does not check for
    overlapping/intersecting triangles.

    Keeps patches separate by renumbering.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "triSurface.H"
#include "OFstream.H"
#include "IFstream.H"
#include "triFace.H"
#include "triFaceList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    argList::addNote
    (
        "add two surfaces via a geometric merge on points."
    );

    argList::validArgs.append("surface file");
    argList::validArgs.append("surface file");
    argList::validArgs.append("output surface file");

    argList::addOption
    (
        "points",
        "file",
        "provide additional points"
    );
    argList::addBoolOption
    (
        "mergeRegions",
        "combine regions from both surfaces"
    );

    argList args(argc, argv);

    const fileName inFileName1 = args[1];
    const fileName inFileName2 = args[2];
    const fileName outFileName = args[3];

    const bool addPoint     = args.optionFound("points");
    const bool mergeRegions = args.optionFound("mergeRegions");

    if (addPoint)
    {
        Info<< "Reading a surface and adding points from a file"
            << "; merging the points and writing the surface to another file"
            << nl << endl;

        Info<< "Surface  : " << inFileName1<< nl
            << "Points   : " << args["points"] << nl
            << "Writing  : " << outFileName << nl << endl;
    }
    else
    {
        Info<< "Reading two surfaces"
            << "; merging points and writing the surface to another file"
            << nl << endl;

        if (mergeRegions)
        {
            Info<< "Regions from the two files will get merged" << nl
                << "Do not use this option if you want to keep the regions"
                << " separate" << nl << endl;
        }
        else
        {
            Info<< "Regions from the two files will not get merged" << nl
                << "Regions from " << inFileName2 << " will get offset so"
                << " as not to overlap with the regions in " << inFileName1
                << nl << endl;
        }


        Info<< "Surface1 : " << inFileName1<< nl
            << "Surface2 : " << inFileName2<< nl
            << "Writing  : " << outFileName << nl << endl;
    }

    const triSurface surface1(inFileName1);

    Info<< "Surface1:" << endl;
    surface1.writeStats(Info);
    Info<< endl;

    const pointField& points1 = surface1.points();

    // Final surface
    triSurface combinedSurf;

    if (addPoint)
    {
        IFstream pointsFile(args["points"]);
        pointField extraPoints(pointsFile);

        Info<< "Additional Points:" << extraPoints.size() << endl;

        vectorField pointsAll(points1);
        label pointi = pointsAll.size();
        pointsAll.setSize(pointsAll.size() + extraPoints.size());

        forAll(extraPoints, i)
        {
            pointsAll[pointi++] = extraPoints[i];
        }

        combinedSurf = triSurface(surface1, surface1.patches(), pointsAll);
    }
    else
    {
        const triSurface surface2(inFileName2);

        Info<< "Surface2:" << endl;
        surface2.writeStats(Info);
        Info<< endl;


        // Make new storage
        List<labelledTri> facesAll(surface1.size() + surface2.size());

        const pointField& points2 = surface2.points();

        vectorField pointsAll(points1.size() + points2.size());


        label pointi = 0;
        // Copy points1 into pointsAll
        forAll(points1, point1i)
        {
            pointsAll[pointi++] = points1[point1i];
        }
        // Add surface2 points
        forAll(points2, point2i)
        {
            pointsAll[pointi++] = points2[point2i];
        }


        label trianglei = 0;

        // Copy triangles1 into trianglesAll
        forAll(surface1, facei)
        {
            facesAll[trianglei++] = surface1[facei];
        }
        label nRegions1 = surface1.patches().size();


        if (!mergeRegions)
        {
            Info<< "Surface " << inFileName1 << " has " << nRegions1
                << " regions"
                << nl
                << "All region numbers in " << inFileName2 << " will be offset"
                << " by this amount" << nl << endl;
        }

        // Add (renumbered) surface2 triangles
        forAll(surface2, facei)
        {
            const labelledTri& tri = surface2[facei];

            labelledTri& destTri = facesAll[trianglei++];
            destTri[0] = tri[0] + points1.size();
            destTri[1] = tri[1] + points1.size();
            destTri[2] = tri[2] + points1.size();
            if (mergeRegions)
            {
                destTri.region() = tri.region();
            }
            else
            {
                destTri.region() = tri.region() + nRegions1;
            }
        }

        label nRegions2 = surface2.patches().size();

        geometricSurfacePatchList newPatches;

        if (mergeRegions)
        {
            // Overwrite
            newPatches.setSize(max(nRegions1, nRegions2));

            forAll(surface1.patches(), patchi)
            {
                newPatches[patchi] = surface1.patches()[patchi];
            }
            forAll(surface2.patches(), patchi)
            {
                newPatches[patchi] = surface2.patches()[patchi];
            }
        }
        else
        {
            Info<< "Regions from " << inFileName2 << " have been renumbered:"
                << nl
                << "    old\tnew" << nl;

            for (label regionI = 0; regionI < nRegions2; regionI++)
            {
                Info<< "    " << regionI << '\t' << regionI+nRegions1
                    << nl;
            }
            Info<< nl;

            newPatches.setSize(nRegions1 + nRegions2);

            label newPatchi = 0;

            forAll(surface1.patches(), patchi)
            {
                newPatches[newPatchi++] = surface1.patches()[patchi];
            }

            forAll(surface2.patches(), patchi)
            {
                newPatches[newPatchi++] = surface2.patches()[patchi];
            }
        }


        Info<< "New patches:" << nl;
        forAll(newPatches, patchi)
        {
            Info<< "    " << patchi << '\t' << newPatches[patchi].name() << nl;
        }
        Info<< endl;


        // Construct new surface mesh
        combinedSurf = triSurface(facesAll, newPatches, pointsAll);
    }

    // Merge all common points and do some checks
    combinedSurf.cleanup(true);

    Info<< "Merged surface:" << endl;

    combinedSurf.writeStats(Info);

    Info<< endl;

    Info<< "Writing : " << outFileName << endl;

    // No need to 'group' while writing since all in correct order anyway.
    combinedSurf.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
