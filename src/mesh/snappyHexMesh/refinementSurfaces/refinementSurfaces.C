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

\*---------------------------------------------------------------------------*/

#include "refinementSurfaces.H"
#include "Time.H"
#include "searchableSurfaces.H"
#include "shellSurfaces.H"
#include "triSurfaceMesh.H"
#include "labelPair.H"
#include "searchableSurfacesQueries.H"
#include "UPtrList.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementSurfaces::refinementSurfaces
(
    const searchableSurfaces& allGeometry,
    const dictionary& surfacesDict,
    const label gapLevelIncrement
)
:
    allGeometry_(allGeometry),
    surfaces_(surfacesDict.size()),
    names_(surfacesDict.size()),
    surfZones_(surfacesDict.size()),
    regionOffset_(surfacesDict.size())
{
    // Wildcard specification : loop over all surface, all regions
    // and try to find a match.

    // Count number of surfaces.
    label surfI = 0;
    forAll(allGeometry_.names(), geomI)
    {
        const word& geomName = allGeometry_.names()[geomI];

        if (surfacesDict.found(geomName))
        {
            surfI++;
        }
    }

    // Size lists
    surfaces_.setSize(surfI);
    names_.setSize(surfI);
    surfZones_.setSize(surfI);
    regionOffset_.setSize(surfI);

    labelList globalMinLevel(surfI, 0);
    labelList globalMaxLevel(surfI, 0);
    labelList globalLevelIncr(surfI, 0);
    scalarField globalAngle(surfI, -great);
    PtrList<dictionary> globalPatchInfo(surfI);
    List<Map<label>> regionMinLevel(surfI);
    List<Map<label>> regionMaxLevel(surfI);
    List<Map<label>> regionLevelIncr(surfI);
    List<Map<scalar>> regionAngle(surfI);
    List<Map<autoPtr<dictionary>>> regionPatchInfo(surfI);


    HashSet<word> unmatchedKeys(surfacesDict.toc());

    surfI = 0;
    forAll(allGeometry_.names(), geomI)
    {
        const word& geomName = allGeometry_.names()[geomI];

        const entry* ePtr = surfacesDict.lookupEntryPtr(geomName, false, true);

        if (ePtr)
        {
            const dictionary& dict = ePtr->dict();
            unmatchedKeys.erase(ePtr->keyword());

            names_[surfI] = geomName;
            surfaces_[surfI] = geomI;

            const labelPair refLevel(dict.lookup("level"));
            globalMinLevel[surfI] = refLevel[0];
            globalMaxLevel[surfI] = refLevel[1];
            globalLevelIncr[surfI] = dict.lookupOrDefault
            (
                "gapLevelIncrement",
                gapLevelIncrement
            );

            if
            (
                globalMinLevel[surfI] < 0
             || globalMaxLevel[surfI] < globalMinLevel[surfI]
             || globalLevelIncr[surfI] < 0
            )
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "Illegal level specification for surface "
                    << names_[surfI]
                    << " : minLevel:" << globalMinLevel[surfI]
                    << " maxLevel:" << globalMaxLevel[surfI]
                    << " levelIncrement:" << globalLevelIncr[surfI]
                    << exit(FatalIOError);
            }

            const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

            // Surface zones
            surfZones_.set(surfI, new surfaceZonesInfo(surface, dict));

            // Global perpendicular angle
            if (dict.found("patchInfo"))
            {
                globalPatchInfo.set
                (
                    surfI,
                    dict.subDict("patchInfo").clone()
                );
            }
            dict.readIfPresent("perpendicularAngle", globalAngle[surfI]);

            if (dict.found("regions"))
            {
                const dictionary& regionsDict = dict.subDict("regions");
                const wordList& regionNames = surface.regions();

                forAll(regionNames, regionI)
                {
                    if (regionsDict.found(regionNames[regionI]))
                    {
                        // Get the dictionary for region
                        const dictionary& regionDict = regionsDict.subDict
                        (
                            regionNames[regionI]
                        );

                        const labelPair refLevel(regionDict.lookup("level"));

                        regionMinLevel[surfI].insert(regionI, refLevel[0]);
                        regionMaxLevel[surfI].insert(regionI, refLevel[1]);
                        label levelIncr = regionDict.lookupOrDefault
                        (
                            "gapLevelIncrement",
                            gapLevelIncrement
                        );
                        regionLevelIncr[surfI].insert(regionI, levelIncr);

                        if
                        (
                            refLevel[0] < 0
                         || refLevel[1] < refLevel[0]
                         || levelIncr < 0
                        )
                        {
                            FatalIOErrorInFunction
                            (
                                dict
                            )   << "Illegal level specification for surface "
                                << names_[surfI] << " region "
                                << regionNames[regionI]
                                << " : minLevel:" << refLevel[0]
                                << " maxLevel:" << refLevel[1]
                                << " levelIncrement:" << levelIncr
                                << exit(FatalIOError);
                        }

                        if (regionDict.found("perpendicularAngle"))
                        {
                            regionAngle[surfI].insert
                            (
                                regionI,
                                readScalar
                                (
                                    regionDict.lookup("perpendicularAngle")
                                )
                            );
                        }

                        if (regionDict.found("patchInfo"))
                        {
                            regionPatchInfo[surfI].insert
                            (
                                regionI,
                                regionDict.subDict("patchInfo").clone()
                            );
                        }
                    }
                }
            }
            surfI++;
        }
    }

    if (unmatchedKeys.size() > 0)
    {
        IOWarningInFunction
        (
            surfacesDict
        )   << "Not all entries in refinementSurfaces dictionary were used."
            << " The following entries were not used : "
            << unmatchedKeys.sortedToc()
            << endl;
    }


    // Calculate local to global region offset
    label nRegions = 0;

    forAll(surfaces_, surfI)
    {
        regionOffset_[surfI] = nRegions;
        nRegions += allGeometry_[surfaces_[surfI]].regions().size();
    }

    // Rework surface specific information into information per global region
    minLevel_.setSize(nRegions);
    minLevel_ = 0;
    maxLevel_.setSize(nRegions);
    maxLevel_ = 0;
    gapLevel_.setSize(nRegions);
    gapLevel_ = -1;
    perpendicularAngle_.setSize(nRegions);
    perpendicularAngle_ = -great;
    patchInfo_.setSize(nRegions);


    forAll(globalMinLevel, surfI)
    {
        label nRegions = allGeometry_[surfaces_[surfI]].regions().size();

        // Initialise to global (i.e. per surface)
        for (label i = 0; i < nRegions; i++)
        {
            label globalRegionI = regionOffset_[surfI] + i;
            minLevel_[globalRegionI] = globalMinLevel[surfI];
            maxLevel_[globalRegionI] = globalMaxLevel[surfI];
            gapLevel_[globalRegionI] =
                maxLevel_[globalRegionI]
              + globalLevelIncr[surfI];

            perpendicularAngle_[globalRegionI] = globalAngle[surfI];
            if (globalPatchInfo.set(surfI))
            {
                patchInfo_.set
                (
                    globalRegionI,
                    globalPatchInfo[surfI].clone()
                );
            }
        }

        // Overwrite with region specific information
        forAllConstIter(Map<label>, regionMinLevel[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            minLevel_[globalRegionI] = iter();
            maxLevel_[globalRegionI] = regionMaxLevel[surfI][iter.key()];
            gapLevel_[globalRegionI] =
                maxLevel_[globalRegionI]
              + regionLevelIncr[surfI][iter.key()];
        }
        forAllConstIter(Map<scalar>, regionAngle[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            perpendicularAngle_[globalRegionI] = regionAngle[surfI][iter.key()];
        }

        const Map<autoPtr<dictionary>>& localInfo = regionPatchInfo[surfI];
        forAllConstIter(Map<autoPtr<dictionary>>, localInfo, iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            patchInfo_.set(globalRegionI, iter()().clone());
        }
    }
}


Foam::refinementSurfaces::refinementSurfaces
(
    const searchableSurfaces& allGeometry,
    const labelList& surfaces,
    const wordList& names,
    const PtrList<surfaceZonesInfo>& surfZones,
    const labelList& regionOffset,
    const labelList& minLevel,
    const labelList& maxLevel,
    const labelList& gapLevel,
    const scalarField& perpendicularAngle,
    PtrList<dictionary>& patchInfo
)
:
    allGeometry_(allGeometry),
    surfaces_(surfaces),
    names_(names),
    surfZones_(surfZones),
    regionOffset_(regionOffset),
    minLevel_(minLevel),
    maxLevel_(maxLevel),
    gapLevel_(gapLevel),
    perpendicularAngle_(perpendicularAngle),
    patchInfo_(patchInfo.size())
{
    forAll(patchInfo_, pI)
    {
        if (patchInfo.set(pI))
        {
            patchInfo_.set(pI, patchInfo.set(pI, nullptr));
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// // Count number of triangles per surface region
// Foam::labelList Foam::refinementSurfaces::countRegions(const triSurface& s)
// {
//     const geometricSurfacePatchList& regions = s.patches();
//
//     labelList nTris(regions.size(), 0);
//
//     forAll(s, triI)
//     {
//         nTris[s[triI].region()]++;
//     }
//     return nTris;
// }


// // Pre-calculate the refinement level for every element
// void Foam::refinementSurfaces::wantedRefinementLevel
// (
//     const shellSurfaces& shells,
//     const label surfI,
//     const List<pointIndexHit>& info,    // Indices
//     const pointField& ctrs,             // Representative coordinate
//     labelList& minLevelField
// ) const
// {
//     const searchableSurface& geom = allGeometry_[surfaces_[surfI]];
//
//     // Get per element the region
//     labelList region;
//     geom.getRegion(info, region);
//
//     // Initialise fields to region wise minLevel
//     minLevelField.setSize(ctrs.size());
//     minLevelField = -1;
//
//     forAll(minLevelField, i)
//     {
//         if (info[i].hit())
//         {
//             minLevelField[i] = minLevel(surfI, region[i]);
//         }
//     }
//
//     // Find out if triangle inside shell with higher level
//     // What level does shell want to refine fc to?
//     labelList shellLevel;
//     shells.findHigherLevel(ctrs, minLevelField, shellLevel);
//
//     forAll(minLevelField, i)
//     {
//         minLevelField[i] = max(minLevelField[i], shellLevel[i]);
//     }
// }


// Precalculate the refinement level for every element of the searchable
// surface.
void Foam::refinementSurfaces::setMinLevelFields
(
    const shellSurfaces& shells
)
{
    forAll(surfaces_, surfI)
    {
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // Precalculation only makes sense if there are different regions
        // (so different refinement levels possible) and there are some
        // elements. Possibly should have 'enough' elements to have fine
        // enough resolution but for now just make sure we don't catch e.g.
        // searchableBox (size=6)
        if (geom.regions().size() > 1 && geom.globalSize() > 10)
        {
            // Representative local coordinates and bounding sphere
            pointField ctrs;
            scalarField radiusSqr;
            geom.boundingSpheres(ctrs, radiusSqr);

            labelList minLevelField(ctrs.size(), -1);
            {
                // Get the element index in a roundabout way. Problem is e.g.
                // distributed surface where local indices differ from global
                // ones (needed for getRegion call)
                List<pointIndexHit> info;
                geom.findNearest(ctrs, radiusSqr, info);

                // Get per element the region
                labelList region;
                geom.getRegion(info, region);

                // From the region get the surface-wise refinement level
                forAll(minLevelField, i)
                {
                    if (info[i].hit()) // Note: should not be necessary
                    {
                        minLevelField[i] = minLevel(surfI, region[i]);
                    }
                }
            }

            // Find out if triangle inside shell with higher level
            // What level does shell want to refine fc to?
            labelList shellLevel;
            shells.findHigherLevel(ctrs, minLevelField, shellLevel);

            forAll(minLevelField, i)
            {
                minLevelField[i] = max(minLevelField[i], shellLevel[i]);
            }

            // Store minLevelField on surface
            const_cast<searchableSurface&>(geom).setField(minLevelField);
        }
    }
}


// Find intersections of edge. Return -1 or first surface with higher minLevel
// number.
void Foam::refinementSurfaces::findHigherIntersection
(
    const pointField& start,
    const pointField& end,
    const labelList& currentLevel,   // current cell refinement level

    labelList& surfaces,
    labelList& surfaceLevel
) const
{
    surfaces.setSize(start.size());
    surfaces = -1;
    surfaceLevel.setSize(start.size());
    surfaceLevel = -1;

    if (surfaces_.empty())
    {
        return;
    }

    if (surfaces_.size() == 1)
    {
        // Optimisation: single segmented surface. No need to duplicate
        // point storage.

        label surfI = 0;

        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // Do intersection test
        List<pointIndexHit> intersectionInfo(start.size());
        geom.findLineAny(start, end, intersectionInfo);

        // See if a cached level field available
        labelList minLevelField;
        geom.getField(intersectionInfo, minLevelField);
        bool haveLevelField =
        (
            returnReduce(minLevelField.size(), sumOp<label>())
          > 0
        );

        if (!haveLevelField && geom.regions().size() == 1)
        {
            minLevelField = labelList
            (
                intersectionInfo.size(),
                minLevel(surfI, 0)
            );
            haveLevelField = true;
        }

        if (haveLevelField)
        {
            forAll(intersectionInfo, i)
            {
                if
                (
                    intersectionInfo[i].hit()
                 && minLevelField[i] > currentLevel[i]
                )
                {
                    surfaces[i] = surfI;    // index of surface
                    surfaceLevel[i] = minLevelField[i];
                }
            }
            return;
        }
    }



    // Work arrays
    pointField p0(start);
    pointField p1(end);
    labelList intersectionToPoint(identity(start.size()));
    List<pointIndexHit> intersectionInfo(start.size());

    forAll(surfaces_, surfI)
    {
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // Do intersection test
        geom.findLineAny(p0, p1, intersectionInfo);

        // See if a cached level field available
        labelList minLevelField;
        geom.getField(intersectionInfo, minLevelField);

        // Copy all hits into arguments, In-place compact misses.
        label missI = 0;
        forAll(intersectionInfo, i)
        {
            // Get the minLevel for the point
            label minLocalLevel = -1;

            if (intersectionInfo[i].hit())
            {
                // Check if minLevelField for this surface.
                if (minLevelField.size())
                {
                    minLocalLevel = minLevelField[i];
                }
                else
                {
                    // Use the min level for the surface instead. Assume
                    // single region 0.
                    minLocalLevel = minLevel(surfI, 0);
                }
            }


            label pointi = intersectionToPoint[i];

            if (minLocalLevel > currentLevel[pointi])
            {
                // Mark point for refinement
                surfaces[pointi] = surfI;
                surfaceLevel[pointi] = minLocalLevel;
            }
            else
            {
                p0[missI] = start[pointi];
                p1[missI] = end[pointi];
                intersectionToPoint[missI] = pointi;
                missI++;
            }
        }

        // All done? Note that this decision should be synchronised
        if (returnReduce(missI, sumOp<label>()) == 0)
        {
            break;
        }

        // Trim misses
        p0.setSize(missI);
        p1.setSize(missI);
        intersectionToPoint.setSize(missI);
        intersectionInfo.setSize(missI);
    }
}


void Foam::refinementSurfaces::findAllHigherIntersections
(
    const pointField& start,
    const pointField& end,
    const labelList& currentLevel,   // current cell refinement level

    const labelList& globalRegionLevel,

    List<vectorList>& surfaceNormal,
    labelListList& surfaceLevel
) const
{
    surfaceLevel.setSize(start.size());
    surfaceNormal.setSize(start.size());

    if (surfaces_.empty())
    {
        return;
    }

    // Work arrays
    List<List<pointIndexHit>> hitInfo;
    labelList pRegions;
    vectorField pNormals;

    forAll(surfaces_, surfI)
    {
        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        surface.findLineAll(start, end, hitInfo);

        // Repack hits for surface into flat list
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // To avoid overhead of calling getRegion for every point

        label n = 0;
        forAll(hitInfo, pointi)
        {
            n += hitInfo[pointi].size();
        }

        List<pointIndexHit> surfInfo(n);
        labelList pointMap(n);
        n = 0;

        forAll(hitInfo, pointi)
        {
            const List<pointIndexHit>& pHits = hitInfo[pointi];

            forAll(pHits, i)
            {
                surfInfo[n] = pHits[i];
                pointMap[n] = pointi;
                n++;
            }
        }

        labelList surfRegion(n);
        vectorField surfNormal(n);
        surface.getRegion(surfInfo, surfRegion);
        surface.getNormal(surfInfo, surfNormal);

        surfInfo.clear();


        // Extract back into pointwise
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        forAll(surfRegion, i)
        {
            label region = globalRegion(surfI, surfRegion[i]);
            label pointi = pointMap[i];

            if (globalRegionLevel[region] > currentLevel[pointi])
            {
                // Append to pointi info
                label sz = surfaceNormal[pointi].size();
                surfaceNormal[pointi].setSize(sz+1);
                surfaceNormal[pointi][sz] = surfNormal[i];

                surfaceLevel[pointi].setSize(sz+1);
                surfaceLevel[pointi][sz] = globalRegionLevel[region];
            }
        }
    }
}


void Foam::refinementSurfaces::findAllHigherIntersections
(
    const pointField& start,
    const pointField& end,
    const labelList& currentLevel,   // current cell refinement level

    const labelList& globalRegionLevel,

    List<pointList>& surfaceLocation,
    List<vectorList>& surfaceNormal,
    labelListList& surfaceLevel
) const
{
    surfaceLevel.setSize(start.size());
    surfaceNormal.setSize(start.size());
    surfaceLocation.setSize(start.size());

    if (surfaces_.empty())
    {
        return;
    }

    // Work arrays
    List<List<pointIndexHit>> hitInfo;
    labelList pRegions;
    vectorField pNormals;

    forAll(surfaces_, surfI)
    {
        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        surface.findLineAll(start, end, hitInfo);

        // Repack hits for surface into flat list
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // To avoid overhead of calling getRegion for every point

        label n = 0;
        forAll(hitInfo, pointi)
        {
            n += hitInfo[pointi].size();
        }

        List<pointIndexHit> surfInfo(n);
        labelList pointMap(n);
        n = 0;

        forAll(hitInfo, pointi)
        {
            const List<pointIndexHit>& pHits = hitInfo[pointi];

            forAll(pHits, i)
            {
                surfInfo[n] = pHits[i];
                pointMap[n] = pointi;
                n++;
            }
        }

        labelList surfRegion(n);
        vectorField surfNormal(n);
        surface.getRegion(surfInfo, surfRegion);
        surface.getNormal(surfInfo, surfNormal);

        // Extract back into pointwise
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        forAll(surfRegion, i)
        {
            label region = globalRegion(surfI, surfRegion[i]);
            label pointi = pointMap[i];

            if (globalRegionLevel[region] > currentLevel[pointi])
            {
                // Append to pointi info
                label sz = surfaceNormal[pointi].size();
                surfaceLocation[pointi].setSize(sz+1);
                surfaceLocation[pointi][sz] = surfInfo[i].hitPoint();

                surfaceNormal[pointi].setSize(sz+1);
                surfaceNormal[pointi][sz] = surfNormal[i];

                surfaceLevel[pointi].setSize(sz+1);
                surfaceLevel[pointi][sz] = globalRegionLevel[region];
            }
        }
    }
}


void Foam::refinementSurfaces::findNearestIntersection
(
    const labelList& surfacesToTest,
    const pointField& start,
    const pointField& end,

    labelList& surface1,
    List<pointIndexHit>& hit1,
    labelList& region1,
    labelList& surface2,
    List<pointIndexHit>& hit2,
    labelList& region2
) const
{
    // 1. intersection from start to end
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Initialize arguments
    surface1.setSize(start.size());
    surface1 = -1;
    hit1.setSize(start.size());
    region1.setSize(start.size());

    // Current end of segment to test.
    pointField nearest(end);
    // Work array
    List<pointIndexHit> nearestInfo(start.size());
    labelList region;

    forAll(surfacesToTest, testI)
    {
        label surfI = surfacesToTest[testI];

        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        // See if any intersection between start and current nearest
        surface.findLine
        (
            start,
            nearest,
            nearestInfo
        );
        surface.getRegion
        (
            nearestInfo,
            region
        );

        forAll(nearestInfo, pointi)
        {
            if (nearestInfo[pointi].hit())
            {
                hit1[pointi] = nearestInfo[pointi];
                surface1[pointi] = surfI;
                region1[pointi] = region[pointi];
                nearest[pointi] = hit1[pointi].hitPoint();
            }
        }
    }


    // 2. intersection from end to last intersection
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the nearest intersection from end to start. Note that we initialize
    // to the first intersection (if any).
    surface2 = surface1;
    hit2 = hit1;
    region2 = region1;

    // Set current end of segment to test.
    forAll(nearest, pointi)
    {
        if (hit1[pointi].hit())
        {
            nearest[pointi] = hit1[pointi].hitPoint();
        }
        else
        {
            // Disable testing by setting to end.
            nearest[pointi] = end[pointi];
        }
    }

    forAll(surfacesToTest, testI)
    {
        label surfI = surfacesToTest[testI];

        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        // See if any intersection between end and current nearest
        surface.findLine
        (
            end,
            nearest,
            nearestInfo
        );
        surface.getRegion
        (
            nearestInfo,
            region
        );

        forAll(nearestInfo, pointi)
        {
            if (nearestInfo[pointi].hit())
            {
                hit2[pointi] = nearestInfo[pointi];
                surface2[pointi] = surfI;
                region2[pointi] = region[pointi];
                nearest[pointi] = hit2[pointi].hitPoint();
            }
        }
    }


    // Make sure that if hit1 has hit something, hit2 will have at least the
    // same point (due to tolerances it might miss its end point)
    forAll(hit1, pointi)
    {
        if (hit1[pointi].hit() && !hit2[pointi].hit())
        {
            hit2[pointi] = hit1[pointi];
            surface2[pointi] = surface1[pointi];
            region2[pointi] = region1[pointi];
        }
    }
}


void Foam::refinementSurfaces::findNearestIntersection
(
    const labelList& surfacesToTest,
    const pointField& start,
    const pointField& end,

    labelList& surface1,
    List<pointIndexHit>& hit1,
    labelList& region1,
    vectorField& normal1,

    labelList& surface2,
    List<pointIndexHit>& hit2,
    labelList& region2,
    vectorField& normal2
) const
{
    // 1. intersection from start to end
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Initialize arguments
    surface1.setSize(start.size());
    surface1 = -1;
    hit1.setSize(start.size());
    region1.setSize(start.size());
    region1 = -1;
    normal1.setSize(start.size());
    normal1 = Zero;

    // Current end of segment to test.
    pointField nearest(end);
    // Work array
    List<pointIndexHit> nearestInfo(start.size());
    labelList region;
    vectorField normal;

    forAll(surfacesToTest, testI)
    {
        label surfI = surfacesToTest[testI];
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // See if any intersection between start and current nearest
        geom.findLine(start, nearest, nearestInfo);
        geom.getRegion(nearestInfo, region);
        geom.getNormal(nearestInfo, normal);

        forAll(nearestInfo, pointi)
        {
            if (nearestInfo[pointi].hit())
            {
                hit1[pointi] = nearestInfo[pointi];
                surface1[pointi] = surfI;
                region1[pointi] = region[pointi];
                normal1[pointi] = normal[pointi];
                nearest[pointi] = hit1[pointi].hitPoint();
            }
        }
    }


    // 2. intersection from end to last intersection
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the nearest intersection from end to start. Note that we initialize
    // to the first intersection (if any).
    surface2 = surface1;
    hit2 = hit1;
    region2 = region1;
    normal2 = normal1;

    // Set current end of segment to test.
    forAll(nearest, pointi)
    {
        if (hit1[pointi].hit())
        {
            nearest[pointi] = hit1[pointi].hitPoint();
        }
        else
        {
            // Disable testing by setting to end.
            nearest[pointi] = end[pointi];
        }
    }

    forAll(surfacesToTest, testI)
    {
        label surfI = surfacesToTest[testI];
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // See if any intersection between end and current nearest
        geom.findLine(end, nearest, nearestInfo);
        geom.getRegion(nearestInfo, region);
        geom.getNormal(nearestInfo, normal);

        forAll(nearestInfo, pointi)
        {
            if (nearestInfo[pointi].hit())
            {
                hit2[pointi] = nearestInfo[pointi];
                surface2[pointi] = surfI;
                region2[pointi] = region[pointi];
                normal2[pointi] = normal[pointi];
                nearest[pointi] = hit2[pointi].hitPoint();
            }
        }
    }


    // Make sure that if hit1 has hit something, hit2 will have at least the
    // same point (due to tolerances it might miss its end point)
    forAll(hit1, pointi)
    {
        if (hit1[pointi].hit() && !hit2[pointi].hit())
        {
            hit2[pointi] = hit1[pointi];
            surface2[pointi] = surface1[pointi];
            region2[pointi] = region1[pointi];
            normal2[pointi] = normal1[pointi];
        }
    }
}


void Foam::refinementSurfaces::findAnyIntersection
(
    const pointField& start,
    const pointField& end,

    labelList& hitSurface,
    List<pointIndexHit>& hitInfo
) const
{
    searchableSurfacesQueries::findAnyIntersection
    (
        allGeometry_,
        surfaces_,
        start,
        end,
        hitSurface,
        hitInfo
    );
}


void Foam::refinementSurfaces::findNearest
(
    const labelList& surfacesToTest,
    const pointField& samples,
    const  scalarField& nearestDistSqr,
    labelList& hitSurface,
    List<pointIndexHit>& hitInfo
) const
{
    labelList geometries(UIndirectList<label>(surfaces_, surfacesToTest));

    // Do the tests. Note that findNearest returns index in geometries.
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        geometries,
        samples,
        nearestDistSqr,
        hitSurface,
        hitInfo
    );

    // Rework the hitSurface to be surface (i.e. index into surfaces_)
    forAll(hitSurface, pointi)
    {
        if (hitSurface[pointi] != -1)
        {
            hitSurface[pointi] = surfacesToTest[hitSurface[pointi]];
        }
    }
}


void Foam::refinementSurfaces::findNearestRegion
(
    const labelList& surfacesToTest,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& hitSurface,
    labelList& hitRegion
) const
{
    labelList geometries(UIndirectList<label>(surfaces_, surfacesToTest));

    // Do the tests. Note that findNearest returns index in geometries.
    List<pointIndexHit> hitInfo;
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        geometries,
        samples,
        nearestDistSqr,
        hitSurface,
        hitInfo
    );

    // Rework the hitSurface to be surface (i.e. index into surfaces_)
    forAll(hitSurface, pointi)
    {
        if (hitSurface[pointi] != -1)
        {
            hitSurface[pointi] = surfacesToTest[hitSurface[pointi]];
        }
    }

    // Collect the region
    hitRegion.setSize(hitSurface.size());
    hitRegion = -1;

    forAll(surfacesToTest, i)
    {
        label surfI = surfacesToTest[i];

        // Collect hits for surfI
        const labelList localIndices(findIndices(hitSurface, surfI));

        List<pointIndexHit> localHits
        (
            UIndirectList<pointIndexHit>
            (
                hitInfo,
                localIndices
            )
        );

        labelList localRegion;
        allGeometry_[surfaces_[surfI]].getRegion(localHits, localRegion);

        forAll(localIndices, i)
        {
            hitRegion[localIndices[i]] = localRegion[i];
        }
    }
}


void Foam::refinementSurfaces::findNearestRegion
(
    const labelList& surfacesToTest,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& hitSurface,
    List<pointIndexHit>& hitInfo,
    labelList& hitRegion,
    vectorField& hitNormal
) const
{
    labelList geometries(UIndirectList<label>(surfaces_, surfacesToTest));

    // Do the tests. Note that findNearest returns index in geometries.
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        geometries,
        samples,
        nearestDistSqr,
        hitSurface,
        hitInfo
    );

    // Rework the hitSurface to be surface (i.e. index into surfaces_)
    forAll(hitSurface, pointi)
    {
        if (hitSurface[pointi] != -1)
        {
            hitSurface[pointi] = surfacesToTest[hitSurface[pointi]];
        }
    }

    // Collect the region
    hitRegion.setSize(hitSurface.size());
    hitRegion = -1;
    hitNormal.setSize(hitSurface.size());
    hitNormal = Zero;

    forAll(surfacesToTest, i)
    {
        label surfI = surfacesToTest[i];

        // Collect hits for surfI
        const labelList localIndices(findIndices(hitSurface, surfI));

        List<pointIndexHit> localHits
        (
            UIndirectList<pointIndexHit>
            (
                hitInfo,
                localIndices
            )
        );

        // Region
        labelList localRegion;
        allGeometry_[surfaces_[surfI]].getRegion(localHits, localRegion);

        forAll(localIndices, i)
        {
            hitRegion[localIndices[i]] = localRegion[i];
        }

        // Normal
        vectorField localNormal;
        allGeometry_[surfaces_[surfI]].getNormal(localHits, localNormal);

        forAll(localIndices, i)
        {
            hitNormal[localIndices[i]] = localNormal[i];
        }
    }
}


//// Find intersection with max of edge. Return -1 or the surface
//// with the highest maxLevel above currentLevel
//Foam::label Foam::refinementSurfaces::findHighestIntersection
//(
//    const point& start,
//    const point& end,
//    const label currentLevel,   // current cell refinement level
//
//    pointIndexHit& maxHit
//) const
//{
//    // surface with highest maxlevel
//    label maxSurface = -1;
//    // maxLevel of maxSurface
//    label maxLevel = currentLevel;
//
//    forAll(*this, surfI)
//    {
//        pointIndexHit hit = operator[](surfI).findLineAny(start, end);
//
//        if (hit.hit())
//        {
//            const triSurface& s = operator[](surfI);
//
//            label region = globalRegion(surfI, s[hit.index()].region());
//
//            if (maxLevel_[region] > maxLevel)
//            {
//                maxSurface = surfI;
//                maxLevel = maxLevel_[region];
//                maxHit = hit;
//            }
//        }
//    }
//
//    if (maxSurface == -1)
//    {
//        // maxLevel unchanged. No interesting surface hit.
//        maxHit.setMiss();
//    }
//
//    return maxSurface;
//}


void Foam::refinementSurfaces::findInside
(
    const labelList& testSurfaces,
    const pointField& pt,
    labelList& insideSurfaces
) const
{
    insideSurfaces.setSize(pt.size());
    insideSurfaces = -1;

    forAll(testSurfaces, i)
    {
        label surfI = testSurfaces[i];

        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        const surfaceZonesInfo::areaSelectionAlgo selectionMethod =
            surfZones_[surfI].zoneInside();

        if
        (
            selectionMethod != surfaceZonesInfo::INSIDE
         && selectionMethod != surfaceZonesInfo::OUTSIDE
        )
        {
            FatalErrorInFunction
                << "Trying to use surface "
                << surface.name()
                << " which has non-geometric inside selection method "
                << surfaceZonesInfo::areaSelectionAlgoNames[selectionMethod]
                << exit(FatalError);
        }

        if (surface.hasVolumeType())
        {
            List<volumeType> volType;
            surface.getVolumeType(pt, volType);

            forAll(volType, pointi)
            {
                if (insideSurfaces[pointi] == -1)
                {
                    if
                    (
                        (
                            volType[pointi] == volumeType::INSIDE
                         && selectionMethod == surfaceZonesInfo::INSIDE
                        )
                     || (
                            volType[pointi] == volumeType::OUTSIDE
                         && selectionMethod == surfaceZonesInfo::OUTSIDE
                        )
                    )
                    {
                        insideSurfaces[pointi] = surfI;
                    }
                }
            }
        }
    }
}


// ************************************************************************* //
