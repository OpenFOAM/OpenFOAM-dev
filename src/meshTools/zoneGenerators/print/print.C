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

\*---------------------------------------------------------------------------*/

#include "print.H"
#include "polyMesh.H"
#include "boundSphere.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(print, 0);
        addToRunTimeSelectionTable(zoneGenerator, print, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::print::print
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    zoneType_
    (
        zoneTypesAllNames.lookupOrDefault
        (
            "zoneType",
            dict,
            zoneTypesAll::all
        )
    ),
    zoneNames_
    (
        dict.found("zone")
      ? new wordList(1, dict.lookup<word>("zone"))
      : dict.found("zones")
      ? new wordList(dict.lookup<wordList>("zones"))
      : nullptr
    ),
    geometry_(dict.lookupOrDefault<bool>("geometry", true))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::print::~print()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::print::generate() const
{
    // Determine which zone types to print
    const bool printPointZones =
        zoneType_ == zoneTypesAll::point || zoneType_ == zoneTypesAll::all;
    const bool printFaceZones =
        zoneType_ == zoneTypesAll::face || zoneType_ == zoneTypesAll::all;
    const bool printCellZones =
        zoneType_ == zoneTypesAll::cell || zoneType_ == zoneTypesAll::all;

    // Get a list of which zones to print. Either user-defined, or all zones.
    wordList allZoneNames;
    if (!zoneNames_.valid())
    {
        wordHashSet allZoneNameSet;
        if (printPointZones) allZoneNameSet.insert(mesh_.pointZones().toc());
        if (printFaceZones) allZoneNameSet.insert(mesh_.faceZones().toc());
        if (printCellZones) allZoneNameSet.insert(mesh_.cellZones().toc());
        allZoneNames = allZoneNameSet.toc();
    }
    const wordList& zoneNames =
        zoneNames_.valid() ? zoneNames_() : allZoneNames;

    // Reference bits of the mesh
    const pointField& points = mesh_.points();
    const faceList& faces = mesh_.faces();
    const scalarField& magFaceAreas = mesh_.magFaceAreas();
    const pointField& faceCentres = mesh_.faceCentres();
    const cellList& cells = mesh_.cells();
    const scalarField& cellVolumes = mesh_.cellVolumes();
    const pointField& cellCentres = mesh_.cellCentres();

    // Keep a set of points for which to calculate geometry. Store both forward
    // and backward addressing so that the maps can be efficiently re-used
    // across multiple zones.
    labelList pointZonePoints(points.size(), -1);
    DynamicList<label> zonePointPoints;

    // Print each zone in turn ...
    forAll(zoneNames, zonei)
    {
        const word& name = zoneNames[zonei];

        // Determine what types of zones with this name we are printing
        const bool printPointZone =
            printPointZones && mesh_.pointZones().found(name);
        const bool printFaceZone =
            printFaceZones && mesh_.faceZones().found(name);
        const bool printCellZone =
            printCellZones && mesh_.cellZones().found(name);
        const bool printMultiple =
            printPointZone + printFaceZone + printCellZone > 1;

        // Preamble
        Info<< "Zone"
            << (printMultiple ? "s" : "")
            << " '" << name << "' " << (printMultiple ? "have" : "has");

        // Initialise integral geometry
        label pointZoneSize = 0;
        point pointZoneAverage = point::zero;
        scalar faceZoneMagArea = 0;
        point faceZoneCentroid = point::zero;
        scalar cellZoneVolume = 0;
        point cellZoneCentroid = point::zero;

        // Clear the point set
        forAll(zonePointPoints, zonePointi)
        {
            pointZonePoints[zonePointPoints[zonePointi]] = -1;
        }
        zonePointPoints.clear();

        if (printPointZone)
        {
            const pointZone& pz = mesh_.pointZones()[name];
            const label nPz = returnReduce(pz.size(), sumOp<label>());

            // Print the number of points
            Info<< ' ' << nPz << " points";

            if (geometry_ && nPz)
            {
                // Generate the global number of points and the average position
                pointZoneSize = nPz;
                forAll(pz, i)
                {
                    pointZoneAverage += points[pz[i]];
                }
                reduce(pointZoneAverage, sumOp<point>());
                pointZoneAverage /= pointZoneSize;

                // Store the indices of the points in this zone
                forAll(pz, pzPointi)
                {
                    const label pointi = pz[pzPointi];
                    if (pointZonePoints[pointi] != -1) continue;
                    pointZonePoints[pointi] = zonePointPoints.size();
                    zonePointPoints.append(pointi);
                }
            }
        }

        if (printFaceZone)
        {
            const faceZone& fz = mesh_.faceZones()[name];
            const label nFz = returnReduce(fz.size(), sumOp<label>());

            // Print the number of faces
            if (printPointZone && printCellZone) Info<< ",";
            if (printPointZone && !printCellZone) Info<< " and";
            Info<< ' ' << nFz << " faces";

            if (geometry_ && nFz)
            {
                // Generate the global face-zone area magnitude and centroid
                forAll(fz, fzFacei)
                {
                    faceZoneMagArea += magFaceAreas[fz[fzFacei]];
                    faceZoneCentroid +=
                        magFaceAreas[fz[fzFacei]]*faceCentres[fz[fzFacei]];
                }
                reduce(faceZoneMagArea, sumOp<scalar>());
                reduce(faceZoneCentroid, sumOp<point>());
                faceZoneCentroid /= faceZoneMagArea;

                // Store the indices of the points in this zone
                forAll(fz, fzFacei)
                {
                    const label facei = fz[fzFacei];
                    forAll(faces[facei], facePointi)
                    {
                        const label pointi = faces[facei][facePointi];
                        if (pointZonePoints[pointi] != -1) continue;
                        pointZonePoints[pointi] = zonePointPoints.size();
                        zonePointPoints.append(pointi);
                    }
                }
            }
        }

        if (printCellZone)
        {
            const cellZone& cz = mesh_.cellZones()[name];
            const label nCz = returnReduce(cz.size(), sumOp<label>());

            // Print the number of cells
            if (printPointZone || printFaceZone) Info<< " and";
            Info<< ' ' << nCz << " cells";

            if (geometry_ && nCz)
            {
                // Generate the global cell-zone volume and centroid
                forAll(cz, czCelli)
                {
                    cellZoneVolume += cellVolumes[cz[czCelli]];
                    cellZoneCentroid +=
                        cellVolumes[cz[czCelli]]*cellCentres[cz[czCelli]];
                }
                reduce(cellZoneVolume, sumOp<scalar>());
                reduce(cellZoneCentroid, sumOp<point>());
                cellZoneCentroid /= cellZoneVolume;

                // Store the indices of the points in this zone
                forAll(cz, czCelli)
                {
                    const label celli = cz[czCelli];
                    forAll(cells[celli], cellFacei)
                    {
                        const label facei = cells[celli][cellFacei];
                        forAll(faces[facei], facePointi)
                        {
                            const label pointi = faces[facei][facePointi];
                            if (pointZonePoints[pointi] != -1) continue;
                            pointZonePoints[pointi] = zonePointPoints.size();
                            zonePointPoints.append(pointi);
                        }
                    }
                }
            }
        }

        Info<< (geometry_ ? ":" : "") << nl;

        // Continue if we're not printing geometry
        if (!geometry_) continue;

        // Print integral geometry
        if (printPointZone && pointZoneSize)
        {
            Info<< "                points average = "
                << pointZoneAverage << nl;
        }
        if (printFaceZone && faceZoneMagArea > 0)
        {
            Info<< "           faces area/centroid = "
                << faceZoneMagArea << '/' << faceZoneCentroid << nl;
        }
        if (printCellZone && cellZoneVolume > 0)
        {
            Info<< "         cells volume/centroid = "
                << cellZoneVolume << '/' << cellZoneCentroid << nl;
        }

        // Don't print bounds if there aren't any points
        if (returnReduce(zonePointPoints.size(), sumOp<label>()) == 0) continue;

        // Print bounding geometry
        const boundBox bb(points, zonePointPoints, true);
        Info<< "             bound box min/max = "
            << bb.min() << '/' << bb.max() << nl;
        const boundSphere bs =
            boundSphere::global(UIndirectList<point>(points, zonePointPoints));
        Info<< "    bound sphere centre/radius = "
            << bs.c() << '/' << bs.r() << nl;
    }

    Info<< endl;

    return zoneSet();
}


// ************************************************************************* //
