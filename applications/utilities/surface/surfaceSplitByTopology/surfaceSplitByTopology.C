/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

Description

    Strips any baffle parts of a surface.  A baffle region is one which is
    reached by walking from an open edge, and stopping when a multiply connected
    edge is reached.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "triSurface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    argList::validOptions.clear();
    argList::validArgs.append("surface file");
    argList::validArgs.append("output surface file");
    argList args(argc, argv);

    fileName surfFileName(args[1]);
    Info<< "Reading surface from " << surfFileName << endl;

    fileName outFileName(args[2]);
    fileName outFileBaseName = outFileName.lessExt();
    word outExtension = outFileName.ext();

    // Load surface
    triSurface surf(surfFileName);

    bool anyZoneRemoved = false;

    label iterationNo = 0;
    label iterationLimit = 10;

    Info<< "Splitting off baffle parts " << endl;

    do
    {
        anyZoneRemoved = false;

        labelList faceZone;

        const labelListList& edFaces = surf.edgeFaces();
        const labelListList& faceEds = surf.faceEdges();

        boolList multipleEdges(edFaces.size(), false);

        forAll(multipleEdges, i)
        {
            if (edFaces[i].size() > 2)
            {
                multipleEdges[i] = true;
            }
        }

        label nZones = surf.markZones(multipleEdges, faceZone);

        if (nZones < 2)
        {
            break;
        }

        boolList nonBaffle(faceZone.size(), true);
        boolList baffle(faceZone.size(), true);
        labelList pointMap;
        labelList faceMap;


        for (label z = 0; z < nZones; z++)
        {
            bool keepZone = true;

            forAll(faceZone, f)
            {
                if (faceZone[f] == z)
                {
                    forAll(faceEds[f], fe)
                    {
                        if (edFaces[faceEds[f][fe]].size() < 2)
                        {
                            keepZone = false;

                            anyZoneRemoved = true;

                            break;
                        }
                    }
                }

                if (!keepZone)
                {
                    break;
                }
            }

            forAll(faceZone, f)
            {
                if (faceZone[f] == z)
                {
                    nonBaffle[f] = keepZone;
                    baffle[f] = !keepZone;
                }
            }
        }

        Info<< "    Iteration " << iterationNo << endl;

        triSurface baffleSurf = surf.subsetMesh(baffle, pointMap, faceMap);

        if (baffleSurf.size())
        {
            fileName bafflePartFileName =
            outFileBaseName
          + "_bafflePart_"
          + name(iterationNo)
          + "." + outExtension;

            Info<< "    Writing baffle part to " << bafflePartFileName << endl;

            baffleSurf.write(bafflePartFileName);
        }

        surf = surf.subsetMesh(nonBaffle, pointMap, faceMap);

        if (iterationNo == iterationLimit)
        {
            WarningInFunction
            << "Iteration limit of " << iterationLimit << "reached" << endl;
        }

        iterationNo++;

    } while (anyZoneRemoved && iterationNo < iterationLimit);

    Info<< "Writing new surface to " << outFileName << endl;

    surf.write(outFileName);

    labelList faceZone;

    const labelListList& edFaces = surf.edgeFaces();

    boolList multipleEdges(edFaces.size(), false);

    forAll(multipleEdges, i)
    {
        if (edFaces[i].size() > 2)
        {
            multipleEdges[i] = true;
        }
    }

    label nZones = surf.markZones(multipleEdges, faceZone);

    Info<< "Splitting remaining multiply connected parts" << endl;

    for (label z = 0; z < nZones; z++)
    {

        boolList include(faceZone.size(), false);
        labelList pointMap;
        labelList faceMap;

        forAll(faceZone, f)
        {
            if (faceZone[f] == z)
            {
                include[f] = true;
            }
        }

        triSurface zoneSurf = surf.subsetMesh(include, pointMap, faceMap);


        fileName remainingPartFileName =
            outFileBaseName
          + "_multiplePart_"
          + name(z)
          + "." + outExtension;

        Info<< "    Writing multiple part "
            << z << " to " << remainingPartFileName << endl;

        zoneSurf.write(remainingPartFileName);
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
