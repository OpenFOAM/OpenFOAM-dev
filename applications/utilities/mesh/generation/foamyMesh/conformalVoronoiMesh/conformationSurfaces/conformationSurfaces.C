/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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

#include "conformationSurfaces.H"
#include "conformalVoronoiMesh.H"
#include "triSurface.H"
#include "searchableSurfaceFeatures.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(conformationSurfaces, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::conformationSurfaces::hasBoundedVolume
(
    List<volumeType>& referenceVolumeTypes
) const
{
    vector sum(Zero);
    label totalTriangles = 0;

    forAll(surfaces_, s)
    {
        const searchableSurface& surface(allGeometry_[surfaces_[s]]);

        if
        (
            surface.hasVolumeType()
         && (
                normalVolumeTypes_[regionOffset_[s]]
             != extendedFeatureEdgeMesh::BOTH
            )
        )
        {
            pointField pts(1, locationInMesh_);

            List<volumeType> vTypes
            (
                pts.size(),
                volumeType::unknown
            );

            surface.getVolumeType(pts, vTypes);

            referenceVolumeTypes[s] = vTypes[0];

            Info<< "    is "
                << volumeType::names[referenceVolumeTypes[s]]
                << " surface " << surface.name()
                << endl;
        }

        if (isA<triSurface>(surface))
        {
            const triSurface& triSurf = refCast<const triSurface>(surface);

            const pointField& surfPts = triSurf.points();

            Info<< "    Checking " << surface.name() << endl;

            label nBaffles = 0;

            Info<< "        Index = " << surfaces_[s] << endl;
            Info<< "        Offset = " << regionOffset_[s] << endl;

            forAll(triSurf, sI)
            {
                const label patchID =
                    triSurf[sI].region()
                  + regionOffset_[s];

                // Don't include baffle surfaces in the calculation
                if
                (
                    normalVolumeTypes_[patchID]
                 != extendedFeatureEdgeMesh::BOTH
                )
                {
                    sum += triSurf[sI].area(surfPts);
                }
                else
                {
                    nBaffles++;
                }
            }
            Info<< "        has " << nBaffles << " baffles out of "
                << triSurf.size() << " triangles" << endl;

            totalTriangles += triSurf.size();
        }
    }

    Info<< "    Sum of all the surface normals (if near zero, surface is"
        << " probably closed):" << nl
        << "    Note: Does not include baffle surfaces in calculation" << nl
        << "        Sum = " << sum/(totalTriangles + small) << nl
        << "        mag(Sum) = " << mag(sum)/(totalTriangles + small)
        << endl;
}


void Foam::conformationSurfaces::readFeatures
(
    const label surfI,
    const dictionary& featureDict,
    const word& surfaceName,
    label& featureIndex
)
{
    word featureMethod =
        featureDict.lookupOrDefault<word>("featureMethod", "none");

    if (featureMethod == "extendedFeatureEdgeMesh")
    {
        fileName feMeshName(featureDict.lookup("extendedFeatureEdgeMesh"));

        Info<< "    features: " << feMeshName << endl;

        features_.set
        (
            featureIndex,
            new extendedFeatureEdgeMesh
            (
                IOobject
                (
                    feMeshName,
                    runTime_.time().constant(),
                    "extendedFeatureEdgeMesh",
                    runTime_.time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        featureIndex++;
    }
    else if (featureMethod == "extractFeatures")
    {
        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        Info<< "    features: " << surface.name()
            << " of type " << surface.type()
            << ", id: " << featureIndex << endl;

        autoPtr<searchableSurfaceFeatures> ssFeatures
        (
            searchableSurfaceFeatures::New(surface, featureDict)
        );

        if (ssFeatures().hasFeatures())
        {
            features_.set
            (
                featureIndex,
                ssFeatures().features()
            );

            featureIndex++;
        }
        else
        {
            WarningInFunction
                << surface.name() << " of type "
                << surface.type() << " does not have features"
                << endl;
        }
    }
    else if (featureMethod == "none")
    {
        // Currently nothing to do
    }
    else
    {
        FatalErrorInFunction
            << "No valid featureMethod found for surface " << surfaceName
            << nl << "Use \"extendedFeatureEdgeMesh\" "
            << "or \"extractFeatures\"."
            << exit(FatalError);
    }
}

void Foam::conformationSurfaces::readFeatures
(
    const dictionary& featureDict,
    const word& surfaceName,
    label& featureIndex
)
{
    word featureMethod =
        featureDict.lookupOrDefault<word>("featureMethod", "none");

    if (featureMethod == "extendedFeatureEdgeMesh")
    {
        fileName feMeshName(featureDict.lookup("extendedFeatureEdgeMesh"));

        Info<< "    features: " << feMeshName << ", id: " << featureIndex
            << endl;

        features_.set
        (
            featureIndex,
            new extendedFeatureEdgeMesh
            (
                IOobject
                (
                    feMeshName,
                    runTime_.time().constant(),
                    "extendedFeatureEdgeMesh",
                    runTime_.time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        featureIndex++;
    }
    else if (featureMethod == "none")
    {
        // Currently nothing to do
    }
    else
    {
        FatalErrorInFunction
            << "No valid featureMethod found for surface " << surfaceName
            << nl << "Use \"extendedFeatureEdgeMesh\" "
            << "or \"extractFeatures\"."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conformationSurfaces::conformationSurfaces
(
    const Time& runTime,
    Random& rndGen,
    const searchableSurfaces& allGeometry,
    const dictionary& surfaceConformationDict
)
:
    runTime_(runTime),
    allGeometry_(allGeometry),
    features_(),
    locationInMesh_(surfaceConformationDict.lookup("locationInMesh")),
    surfaces_(),
    allGeometryToSurfaces_(),
    normalVolumeTypes_(),
    patchNames_(),
    surfZones_(),
    regionOffset_(),
    patchInfo_(),
    globalBounds_(),
    referenceVolumeTypes_(0)
{
    const dictionary& surfacesDict
    (
        surfaceConformationDict.subDict("geometryToConformTo")
    );

    const dictionary& additionalFeaturesDict
    (
        surfaceConformationDict.subDict("additionalFeatures")
    );


    // Wildcard specification : loop over all surface, all regions
    // and try to find a match.

    // Count number of surfaces.
    label surfI = 0;
    forAll(allGeometry.names(), geomI)
    {
        const word& geomName = allGeometry_.names()[geomI];

        if (surfacesDict.found(geomName))
        {
            surfI++;
        }
    }

    const label nAddFeat = additionalFeaturesDict.size();

    Info<< nl << "Reading geometryToConformTo" << endl;

    allGeometryToSurfaces_.setSize(allGeometry_.size(), -1);

    normalVolumeTypes_.setSize(surfI);
    surfaces_.setSize(surfI);
    surfZones_.setSize(surfI);

    // Features may be attached to host surfaces or independent
    features_.setSize(surfI + nAddFeat);

    label featureI = 0;

    regionOffset_.setSize(surfI, 0);

    PtrList<dictionary> globalPatchInfo(surfI);
    List<Map<autoPtr<dictionary>>> regionPatchInfo(surfI);
    List<sideVolumeType> globalVolumeTypes(surfI);
    List<Map<sideVolumeType>> regionVolumeTypes(surfI);

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

            surfaces_[surfI] = geomI;

            const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

            // Surface zones
            if (dict.found("faceZone"))
            {
                surfZones_.set(surfI, new surfaceZonesInfo(surface, dict));
            }

            allGeometryToSurfaces_[surfaces_[surfI]] = surfI;

            Info<< nl << "    " << geomName << endl;

            const wordList& regionNames =
                allGeometry_.regionNames()[surfaces_[surfI]];

            patchNames_.append(regionNames);

            globalVolumeTypes[surfI] =
            (
                extendedFeatureEdgeMesh::sideVolumeTypeNames_
                [
                    dict.lookupOrDefault<word>
                    (
                        "meshableSide",
                        "inside"
                    )
                ]
            );

            if (!globalVolumeTypes[surfI])
            {
                if (!surface.hasVolumeType())
                {
                    WarningInFunction
                        << "Non-baffle surface "
                        << surface.name()
                        << " does not allow inside/outside queries."
                        << " This usually is an error." << endl;
                }
            }

            // Load patch info
            if (dict.found("patchInfo"))
            {
                globalPatchInfo.set
                (
                    surfI,
                    dict.subDict("patchInfo").clone()
                );
            }

            readFeatures
            (
                surfI,
                dict,
                geomName,
                featureI
            );

            const wordList& rNames = surface.regions();

            if (dict.found("regions"))
            {
                const dictionary& regionsDict = dict.subDict("regions");

                forAll(rNames, regionI)
                {
                    const word& regionName = rNames[regionI];

                    if (regionsDict.found(regionName))
                    {
                        Info<< "        region " << regionName << endl;

                        // Get the dictionary for region
                        const dictionary& regionDict = regionsDict.subDict
                        (
                            regionName
                        );

                        if (regionDict.found("patchInfo"))
                        {
                            regionPatchInfo[surfI].insert
                            (
                                regionI,
                                regionDict.subDict("patchInfo").clone()
                            );
                        }

                        regionVolumeTypes[surfI].insert
                        (
                            regionI,
                            extendedFeatureEdgeMesh::sideVolumeTypeNames_
                            [
                                 regionDict.lookupOrDefault<word>
                                 (
                                     "meshableSide",
                                     extendedFeatureEdgeMesh::
                                     sideVolumeTypeNames_
                                     [
                                        globalVolumeTypes[surfI]
                                     ]
                                 )
                            ]
                        );

                        readFeatures(regionDict, regionName, featureI);
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
        )   << "Not all entries in conformationSurfaces dictionary were used."
            << " The following entries were not used : "
            << unmatchedKeys.sortedToc()
            << endl;
    }


    // Calculate local to global region offset
    label nRegions = 0;

    forAll(surfaces_, surfI)
    {
        regionOffset_[surfI] = nRegions;

        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];
        nRegions += surface.regions().size();
    }

    // Rework surface specific information into information per global region
    patchInfo_.setSize(nRegions);
    normalVolumeTypes_.setSize(nRegions);

    forAll(surfaces_, surfI)
    {
        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        label nRegions = surface.regions().size();

        // Initialise to global (i.e. per surface)
        for (label i = 0; i < nRegions; i++)
        {
            label globalRegionI = regionOffset_[surfI] + i;
            normalVolumeTypes_[globalRegionI] = globalVolumeTypes[surfI];
            if (globalPatchInfo.set(surfI))
            {
                patchInfo_.set
                (
                    globalRegionI,
                    globalPatchInfo[surfI].clone()
                );
            }
        }

        forAllConstIter(Map<sideVolumeType>, regionVolumeTypes[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            normalVolumeTypes_[globalRegionI] =
                regionVolumeTypes[surfI][iter.key()];
        }

        const Map<autoPtr<dictionary>>& localInfo = regionPatchInfo[surfI];
        forAllConstIter(Map<autoPtr<dictionary>>, localInfo, iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            patchInfo_.set(globalRegionI, iter()().clone());
        }
    }



    if (!additionalFeaturesDict.empty())
    {
        Info<< nl << "Reading additionalFeatures" << endl;
    }

    forAllConstIter(dictionary, additionalFeaturesDict, iter)
    {
        word featureName = iter().keyword();

        Info<< nl << "    " << iter().keyword() << endl;

        const dictionary& featureSubDict
        (
            additionalFeaturesDict.subDict(featureName)
        );

        readFeatures(featureSubDict, featureName, featureI);
    }

    // Remove unnecessary space from the features list
    features_.setSize(featureI);

    globalBounds_ = treeBoundBox
    (
        searchableSurfacesQueries::bounds(allGeometry_, surfaces_)
    );

    // Extend the global bounds to stop the bound box sitting on the surfaces
    // to be conformed to
    vector newSpan = 1e-4*globalBounds_.span();
    globalBounds_.min() -= newSpan;
    globalBounds_.max() += newSpan;

    // Look at all surfaces at determine whether the locationInMesh point is
    // inside or outside each, to establish a signature for the domain to be
    // meshed.

    referenceVolumeTypes_.setSize
    (
        surfaces_.size(),
        volumeType::unknown
    );

    Info<< endl
        << "Testing for locationInMesh " << locationInMesh_ << endl;

    hasBoundedVolume(referenceVolumeTypes_);

    if (debug)
    {
        Info<< "Names = " << allGeometry_.names() << endl;
        Info<< "Surfaces = " << surfaces_ << endl;
        Info<< "AllGeom to Surfaces = " << allGeometryToSurfaces_ << endl;
        Info<< "Volume types = " << normalVolumeTypes_ << endl;
        Info<< "Patch names = " << patchNames_ << endl;
        Info<< "Region Offset = " << regionOffset_ << endl;

        forAll(features_, fI)
        {
            Info<< features_[fI].name() << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conformationSurfaces::~conformationSurfaces()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::conformationSurfaces::overlaps(const treeBoundBox& bb) const
{
    forAll(surfaces_, s)
    {
        if (allGeometry_[surfaces_[s]].overlaps(bb))
        {
            return true;
        }
    }

    return false;
}


Foam::Field<bool> Foam::conformationSurfaces::inside
(
    const pointField& samplePts
) const
{
    return wellInside(samplePts, scalarField(samplePts.size(), 0.0));
}


bool Foam::conformationSurfaces::inside
(
    const point& samplePt
) const
{
    return wellInside(pointField(1, samplePt), scalarField(1, 0))[0];
}


Foam::Field<bool> Foam::conformationSurfaces::outside
(
    const pointField& samplePts
) const
{
    return wellOutside(samplePts, scalarField(samplePts.size(), 0.0));
}


bool Foam::conformationSurfaces::outside
(
    const point& samplePt
) const
{
    return wellOutside(pointField(1, samplePt), scalarField(1, 0))[0];
    // return !inside(samplePt);
}


Foam::Field<bool> Foam::conformationSurfaces::wellInOutSide
(
    const pointField& samplePts,
    const scalarField& testDistSqr,
    const bool testForInside
) const
{
    List<List<volumeType>> surfaceVolumeTests
    (
        surfaces_.size(),
        List<volumeType>
        (
            samplePts.size(),
            volumeType::unknown
        )
    );

    // Get lists for the volumeTypes for each sample wrt each surface
    forAll(surfaces_, s)
    {
        const searchableSurface& surface(allGeometry_[surfaces_[s]]);

        const label regionI = regionOffset_[s];

        if (normalVolumeTypes_[regionI] != extendedFeatureEdgeMesh::BOTH)
        {
            surface.getVolumeType(samplePts, surfaceVolumeTests[s]);
        }
    }

    // Compare the volumeType result for each point wrt to each surface with the
    // reference value and if the points are inside the surface by a given
    // distanceSquared

    // Assume that the point is wellInside until demonstrated otherwise.
    Field<bool> insideOutsidePoint(samplePts.size(), testForInside);

    // Check if the points are inside the surface by the given distance squared

    labelList hitSurfaces;
    List<pointIndexHit> hitInfo;
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        surfaces_,
        samplePts,
        testDistSqr,
        hitSurfaces,
        hitInfo
    );

    forAll(samplePts, i)
    {
        const pointIndexHit& pHit = hitInfo[i];

        if (pHit.hit())
        {
            // If the point is within range of the surface, then it can't be
            // well (in|out)side
            insideOutsidePoint[i] = false;

            continue;
        }

        forAll(surfaces_, s)
        {
            const label regionI = regionOffset_[s];

            if (normalVolumeTypes_[regionI] == extendedFeatureEdgeMesh::BOTH)
            {
                continue;
            }

            const searchableSurface& surface(allGeometry_[surfaces_[s]]);

            if
            (
                !surface.hasVolumeType()
             //&& surfaceVolumeTests[s][i] == volumeType::unknown
            )
            {
                pointField sample(1, samplePts[i]);
                scalarField nearestDistSqr(1, great);
                List<pointIndexHit> info;

                surface.findNearest(sample, nearestDistSqr, info);

                vector hitDir = info[0].rawPoint() - samplePts[i];
                hitDir /= mag(hitDir) + small;

                pointIndexHit surfHit;
                label hitSurface;

                findSurfaceNearestIntersection
                (
                    samplePts[i],
                    info[0].rawPoint() - 1e-3*mag(hitDir)*hitDir,
                    surfHit,
                    hitSurface
                );

                if (surfHit.hit() && hitSurface != surfaces_[s])
                {
                    continue;
                }
            }

            if (surfaceVolumeTests[s][i] == volumeType::outside)
            {
                if
                (
                    normalVolumeTypes_[regionI]
                 == extendedFeatureEdgeMesh::INSIDE
                )
                {
                    insideOutsidePoint[i] = !testForInside;
                    break;
                }
            }
            else if (surfaceVolumeTests[s][i] == volumeType::inside)
            {
                if
                (
                    normalVolumeTypes_[regionI]
                 == extendedFeatureEdgeMesh::OUTSIDE
                )
                {
                    insideOutsidePoint[i] = !testForInside;
                    break;
                }
            }
        }
    }

    return insideOutsidePoint;
}


Foam::Field<bool> Foam::conformationSurfaces::wellInside
(
    const pointField& samplePts,
    const scalarField& testDistSqr
) const
{
    return wellInOutSide(samplePts, testDistSqr, true);
}


bool Foam::conformationSurfaces::wellInside
(
    const point& samplePt,
    scalar testDistSqr
) const
{
    return wellInside(pointField(1, samplePt), scalarField(1, testDistSqr))[0];
}


Foam::Field<bool> Foam::conformationSurfaces::wellOutside
(
    const pointField& samplePts,
    const scalarField& testDistSqr
) const
{
    return wellInOutSide(samplePts, testDistSqr, false);
}


bool Foam::conformationSurfaces::wellOutside
(
    const point& samplePt,
    scalar testDistSqr
) const
{
    return wellOutside(pointField(1, samplePt), scalarField(1, testDistSqr))[0];
}


bool Foam::conformationSurfaces::findSurfaceAnyIntersection
(
    const point& start,
    const point& end
) const
{
    labelList hitSurfaces;
    List<pointIndexHit> hitInfo;

    searchableSurfacesQueries::findAnyIntersection
    (
        allGeometry_,
        surfaces_,
        pointField(1, start),
        pointField(1, end),
        hitSurfaces,
        hitInfo
    );

    return hitInfo[0].hit();
}


void Foam::conformationSurfaces::findSurfaceAnyIntersection
(
    const point& start,
    const point& end,
    pointIndexHit& surfHit,
    label& hitSurface
) const
{
    labelList hitSurfaces;
    List<pointIndexHit> hitInfo;

    searchableSurfacesQueries::findAnyIntersection
    (
        allGeometry_,
        surfaces_,
        pointField(1, start),
        pointField(1, end),
        hitSurfaces,
        hitInfo
    );

    surfHit = hitInfo[0];

    if (surfHit.hit())
    {
        // hitSurfaces has returned the index of the entry in surfaces_ that was
        // found, not the index of the surface in allGeometry_, translating this
        // to allGeometry_

        hitSurface = surfaces_[hitSurfaces[0]];
    }
}


void Foam::conformationSurfaces::findSurfaceAllIntersections
(
    const point& start,
    const point& end,
    List<pointIndexHit>& surfHit,
    labelList& hitSurface
) const
{
    labelListList hitSurfaces;
    List<List<pointIndexHit>> hitInfo;

    searchableSurfacesQueries::findAllIntersections
    (
        allGeometry_,
        surfaces_,
        pointField(1, start),
        pointField(1, end),
        hitSurfaces,
        hitInfo
    );

    surfHit = hitInfo[0];

    hitSurface.setSize(hitSurfaces[0].size());

    forAll(hitSurfaces[0], surfI)
    {
        // hitSurfaces has returned the index of the entry in surfaces_ that was
        // found, not the index of the surface in allGeometry_, translating this
        // to allGeometry_

        hitSurface[surfI] = surfaces_[hitSurfaces[0][surfI]];
    }
}


void Foam::conformationSurfaces::findSurfaceNearestIntersection
(
    const point& start,
    const point& end,
    pointIndexHit& surfHit,
    label& hitSurface
) const
{
    labelList hitSurfacesStart;
    List<pointIndexHit> hitInfoStart;
    labelList hitSurfacesEnd;
    List<pointIndexHit> hitInfoEnd;

    searchableSurfacesQueries::findNearestIntersection
    (
        allGeometry_,
        surfaces_,
        pointField(1, start),
        pointField(1, end),
        hitSurfacesStart,
        hitInfoStart,
        hitSurfacesEnd,
        hitInfoEnd
    );

    surfHit = hitInfoStart[0];

    if (surfHit.hit())
    {
        // hitSurfaces has returned the index of the entry in surfaces_ that was
        // found, not the index of the surface in allGeometry_, translating this
        // to allGeometry_

        hitSurface = surfaces_[hitSurfacesStart[0]];
    }
}


void Foam::conformationSurfaces::findSurfaceNearest
(
    const point& sample,
    scalar nearestDistSqr,
    pointIndexHit& surfHit,
    label& hitSurface
) const
{
    labelList hitSurfaces;
    List<pointIndexHit> surfaceHits;

    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        surfaces_,
        pointField(1, sample),
        scalarField(1, nearestDistSqr),
        hitSurfaces,
        surfaceHits
    );

    surfHit = surfaceHits[0];

    if (surfHit.hit())
    {
        // hitSurfaces has returned the index of the entry in surfaces_ that was
        // found, not the index of the surface in allGeometry_, translating this
        // to allGeometry_

        hitSurface = surfaces_[hitSurfaces[0]];
    }
}


void Foam::conformationSurfaces::findSurfaceNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& surfaceHits,
    labelList& hitSurfaces
) const
{
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        surfaces_,
        samples,
        nearestDistSqr,
        hitSurfaces,
        surfaceHits
    );

    forAll(surfaceHits, i)
    {
        if (surfaceHits[i].hit())
        {
            // hitSurfaces has returned the index of the entry in surfaces_ that
            // was found, not the index of the surface in allGeometry_,
            // translating this to the surface in allGeometry_.

            hitSurfaces[i] = surfaces_[hitSurfaces[i]];
        }
    }
}


void Foam::conformationSurfaces::findFeaturePointNearest
(
    const point& sample,
    scalar nearestDistSqr,
    pointIndexHit& fpHit,
    label& featureHit
) const
{
    // Work arrays
    scalar minDistSqr = nearestDistSqr;
    pointIndexHit hitInfo;

    forAll(features_, testI)
    {
        features_[testI].nearestFeaturePoint
        (
            sample,
            minDistSqr,
            hitInfo
        );

        if (hitInfo.hit())
        {
            minDistSqr = magSqr(hitInfo.hitPoint()- sample);
            fpHit = hitInfo;
            featureHit = testI;
        }
    }
}


void Foam::conformationSurfaces::findEdgeNearest
(
    const point& sample,
    scalar nearestDistSqr,
    pointIndexHit& edgeHit,
    label& featureHit
) const
{
    pointField samples(1, sample);
    scalarField nearestDistsSqr(1, nearestDistSqr);

    List<pointIndexHit> edgeHits;
    labelList featuresHit;

    findEdgeNearest
    (
        samples,
        nearestDistsSqr,
        edgeHits,
        featuresHit
    );

    edgeHit = edgeHits[0];
    featureHit = featuresHit[0];
}


void Foam::conformationSurfaces::findEdgeNearest
(
    const pointField& samples,
    const scalarField& nearestDistsSqr,
    List<pointIndexHit>& edgeHits,
    labelList& featuresHit
) const
{
    // Initialise
    featuresHit.setSize(samples.size());
    featuresHit = -1;
    edgeHits.setSize(samples.size());

    // Work arrays
    scalarField minDistSqr(nearestDistsSqr);
    List<pointIndexHit> hitInfo(samples.size());

    forAll(features_, testI)
    {
        features_[testI].nearestFeatureEdge
        (
            samples,
            minDistSqr,
            hitInfo
        );

        // Update minDistSqr and arguments
        forAll(hitInfo, pointi)
        {
            if (hitInfo[pointi].hit())
            {
                minDistSqr[pointi] = magSqr
                (
                    hitInfo[pointi].hitPoint()
                  - samples[pointi]
                );
                edgeHits[pointi] = hitInfo[pointi];
                featuresHit[pointi] = testI;
            }
        }
    }
}


void Foam::conformationSurfaces::findEdgeNearestByType
(
    const point& sample,
    scalar nearestDistSqr,
    List<pointIndexHit>& edgeHits,
    List<label>& featuresHit
) const
{
    // Initialise
    featuresHit.setSize(extendedFeatureEdgeMesh::nEdgeTypes);
    featuresHit = -1;
    edgeHits.setSize(extendedFeatureEdgeMesh::nEdgeTypes);

    // Work arrays
    scalarField minDistSqr(extendedFeatureEdgeMesh::nEdgeTypes, nearestDistSqr);
    List<pointIndexHit> hitInfo(extendedFeatureEdgeMesh::nEdgeTypes);

    forAll(features_, testI)
    {
        features_[testI].nearestFeatureEdgeByType
        (
            sample,
            minDistSqr,
            hitInfo
        );

        // Update minDistSqr and arguments
        forAll(hitInfo, typeI)
        {
            if (hitInfo[typeI].hit())
            {
                minDistSqr[typeI] = magSqr(hitInfo[typeI].hitPoint() - sample);
                edgeHits[typeI] = hitInfo[typeI];
                featuresHit[typeI] = testI;
            }
        }
    }
}


void Foam::conformationSurfaces::findAllNearestEdges
(
    const point& sample,
    const scalar searchRadiusSqr,
    List<List<pointIndexHit>>& edgeHitsByFeature,
    List<label>& featuresHit
) const
{
    // Initialise
    // featuresHit.setSize(features_.size());
    // featuresHit = -1;
    // edgeHitsByFeature.setSize(features_.size());

    // Work arrays
    List<pointIndexHit> hitInfo(extendedFeatureEdgeMesh::nEdgeTypes);

    forAll(features_, testI)
    {
        features_[testI].allNearestFeatureEdges
        (
            sample,
            searchRadiusSqr,
            hitInfo
        );

        bool anyHit = false;
        forAll(hitInfo, hitI)
        {
            if (hitInfo[hitI].hit())
            {
                anyHit = true;
            }
        }

        if (anyHit)
        {
            edgeHitsByFeature.append(hitInfo);
            featuresHit.append(testI);
        }
    }
}


void Foam::conformationSurfaces::writeFeatureObj(const fileName& prefix) const
{
    OFstream ftStr(runTime_.time().path()/prefix + "_allFeatures.obj");

    Pout<< nl << "Writing all features to " << ftStr.name() << endl;

    label verti = 0;

    forAll(features_, i)
    {
        const extendedFeatureEdgeMesh& fEM(features_[i]);
        const pointField pts(fEM.points());
        const edgeList eds(fEM.edges());

        ftStr << "g " << fEM.name() << endl;

        forAll(eds, j)
        {
            const edge& e = eds[j];

            meshTools::writeOBJ(ftStr, pts[e[0]]); verti++;
            meshTools::writeOBJ(ftStr, pts[e[1]]); verti++;
            ftStr << "l " << verti-1 << ' ' << verti << endl;
        }
    }
}


Foam::label Foam::conformationSurfaces::findPatch
(
    const point& ptA,
    const point& ptB
) const
{
    pointIndexHit surfHit;
    label hitSurface;

    findSurfaceAnyIntersection(ptA, ptB, surfHit, hitSurface);

    return getPatchID(hitSurface, surfHit);
}


Foam::label Foam::conformationSurfaces::findPatch(const point& pt) const
{
    pointIndexHit surfHit;
    label hitSurface;

    findSurfaceNearest(pt, sqr(great), surfHit, hitSurface);

    return getPatchID(hitSurface, surfHit);
}


Foam::label Foam::conformationSurfaces::getPatchID
(
    const label hitSurface,
    const pointIndexHit& surfHit
) const
{
    if (!surfHit.hit())
    {
        return -1;
    }

    labelList surfLocalRegion;

    allGeometry_[hitSurface].getRegion
    (
        List<pointIndexHit>(1, surfHit),
        surfLocalRegion
    );

    const label patchID =
        surfLocalRegion[0]
      + regionOffset_[allGeometryToSurfaces_[hitSurface]];

    return patchID;
}


Foam::extendedFeatureEdgeMesh::sideVolumeType
Foam::conformationSurfaces::meshableSide
(
    const label hitSurface,
    const pointIndexHit& surfHit
) const
{
    const label patchID = getPatchID(hitSurface, surfHit);

    if (patchID == -1)
    {
        return extendedFeatureEdgeMesh::NEITHER;
    }

    return normalVolumeTypes_[patchID];
}


void Foam::conformationSurfaces::getNormal
(
    const label hitSurface,
    const List<pointIndexHit>& surfHit,
    vectorField& normal
) const
{
    allGeometry_[hitSurface].getNormal(surfHit, normal);

    const label patchID = regionOffset_[allGeometryToSurfaces_[hitSurface]];

    // Now flip sign of normal depending on mesh side
    if (normalVolumeTypes_[patchID] == extendedFeatureEdgeMesh::OUTSIDE)
    {
        normal *= -1;
    }
}


// ************************************************************************* //
