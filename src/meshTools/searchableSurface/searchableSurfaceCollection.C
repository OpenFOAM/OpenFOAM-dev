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

\*---------------------------------------------------------------------------*/

#include "searchableSurfaceCollection.H"
#include "addToRunTimeSelectionTable.H"
#include "SortableList.H"
#include "Time.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(searchableSurfaceCollection, 0);
addToRunTimeSelectionTable
(
    searchableSurface,
    searchableSurfaceCollection,
    dict
);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::searchableSurfaceCollection::findNearest
(
    const pointField& samples,
    scalarField& minDistSqr,
    List<pointIndexHit>& nearestInfo,
    labelList& nearestSurf
) const
{
    // Initialise
    nearestInfo.setSize(samples.size());
    nearestInfo = pointIndexHit();
    nearestSurf.setSize(samples.size());
    nearestSurf = -1;

    List<pointIndexHit> hitInfo(samples.size());

    const scalarField localMinDistSqr(samples.size(), GREAT);

    forAll(subGeom_, surfI)
    {
        subGeom_[surfI].findNearest
        (
            cmptDivide  // Transform then divide
            (
                transform_[surfI].localPosition(samples),
                scale_[surfI]
            ),
            localMinDistSqr,
            hitInfo
        );

        forAll(hitInfo, pointI)
        {
            if (hitInfo[pointI].hit())
            {
                // Rework back into global coordinate sys. Multiply then
                // transform
                point globalPt = transform_[surfI].globalPosition
                (
                    cmptMultiply
                    (
                        hitInfo[pointI].rawPoint(),
                        scale_[surfI]
                    )
                );

                scalar distSqr = magSqr(globalPt - samples[pointI]);

                if (distSqr < minDistSqr[pointI])
                {
                    minDistSqr[pointI] = distSqr;
                    nearestInfo[pointI].setPoint(globalPt);
                    nearestInfo[pointI].setHit();
                    nearestInfo[pointI].setIndex
                    (
                        hitInfo[pointI].index()
                      + indexOffset_[surfI]
                    );
                    nearestSurf[pointI] = surfI;
                }
            }
        }
    }
}


// Sort hits into per-surface bins. Misses are rejected. Maintains map back
// to position
void Foam::searchableSurfaceCollection::sortHits
(
    const List<pointIndexHit>& info,
    List<List<pointIndexHit> >& surfInfo,
    labelListList& infoMap
) const
{
    // Count hits per surface.
    labelList nHits(subGeom_.size(), 0);

    forAll(info, pointI)
    {
        if (info[pointI].hit())
        {
            label index = info[pointI].index();
            label surfI = findLower(indexOffset_, index+1);
            nHits[surfI]++;
        }
    }

    // Per surface the hit
    surfInfo.setSize(subGeom_.size());
    // Per surface the original position
    infoMap.setSize(subGeom_.size());

    forAll(surfInfo, surfI)
    {
        surfInfo[surfI].setSize(nHits[surfI]);
        infoMap[surfI].setSize(nHits[surfI]);
    }
    nHits = 0;

    forAll(info, pointI)
    {
        if (info[pointI].hit())
        {
            label index = info[pointI].index();
            label surfI = findLower(indexOffset_, index+1);

            // Store for correct surface and adapt indices back to local
            // ones
            label localI = nHits[surfI]++;
            surfInfo[surfI][localI] = pointIndexHit
            (
                info[pointI].hit(),
                info[pointI].rawPoint(),
                index-indexOffset_[surfI]
            );
            infoMap[surfI][localI] = pointI;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceCollection::searchableSurfaceCollection
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    instance_(dict.size()),
    scale_(dict.size()),
    transform_(dict.size()),
    subGeom_(dict.size()),
    mergeSubRegions_(dict.lookup("mergeSubRegions")),
    indexOffset_(dict.size()+1)
{
    Info<< "SearchableCollection : " << name() << endl;

    label surfI = 0;
    label startIndex = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (dict.isDict(iter().keyword()))
        {
            instance_[surfI] = iter().keyword();

            const dictionary& subDict = dict.subDict(instance_[surfI]);

            scale_[surfI] = subDict.lookup("scale");
            transform_.set
            (
                surfI,
                coordinateSystem::New
                (
                    subDict.subDict("transform")
                )
            );

            const word subGeomName(subDict.lookup("surface"));
            //Pout<< "Trying to find " << subGeomName << endl;

            const searchableSurface& s =
                io.db().lookupObject<searchableSurface>(subGeomName);

            // I don't know yet how to handle the globalSize combined with
            // regionOffset. Would cause non-consecutive indices locally
            // if all indices offset by globalSize() of the local region...
            if (s.size() != s.globalSize())
            {
                FatalErrorIn
                (
                    "searchableSurfaceCollection::searchableSurfaceCollection"
                    "(const IOobject&, const dictionary&)"
                )   << "Cannot use a distributed surface in a collection."
                    << exit(FatalError);
            }

            subGeom_.set(surfI, &const_cast<searchableSurface&>(s));

            indexOffset_[surfI] = startIndex;
            startIndex += subGeom_[surfI].size();

            Info<< "    instance : " << instance_[surfI] << endl;
            Info<< "    surface  : " << s.name() << endl;
            Info<< "    scale    : " << scale_[surfI] << endl;
            Info<< "    coordsys : " << transform_[surfI] << endl;

            surfI++;
        }
    }
    indexOffset_[surfI] = startIndex;

    instance_.setSize(surfI);
    scale_.setSize(surfI);
    transform_.setSize(surfI);
    subGeom_.setSize(surfI);
    indexOffset_.setSize(surfI+1);

    // Bounds is the overall bounds
    bounds() = boundBox(point::max, point::min);

    forAll(subGeom_, surfI)
    {
        const boundBox& surfBb = subGeom_[surfI].bounds();

        // Transform back to global coordinate sys.
        const point surfBbMin = transform_[surfI].globalPosition
        (
            cmptMultiply
            (
                surfBb.min(),
                scale_[surfI]
            )
        );
        const point surfBbMax = transform_[surfI].globalPosition
        (
            cmptMultiply
            (
                surfBb.max(),
                scale_[surfI]
            )
        );

        bounds().min() = min(bounds().min(), surfBbMin);
        bounds().max() = max(bounds().max(), surfBbMax);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableSurfaceCollection::~searchableSurfaceCollection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchableSurfaceCollection::regions() const
{
    if (regions_.size() == 0)
    {
        regionOffset_.setSize(subGeom_.size());

        DynamicList<word> allRegions;
        forAll(subGeom_, surfI)
        {
            regionOffset_[surfI] = allRegions.size();

            if (mergeSubRegions_)
            {
                // Single name regardless how many regions subsurface has
                allRegions.append(instance_[surfI] + "_" + Foam::name(surfI));
            }
            else
            {
                const wordList& subRegions = subGeom_[surfI].regions();

                forAll(subRegions, i)
                {
                    allRegions.append(instance_[surfI] + "_" + subRegions[i]);
                }
            }
        }
        regions_.transfer(allRegions.shrink());
    }
    return regions_;
}


Foam::label Foam::searchableSurfaceCollection::size() const
{
    return indexOffset_.last();
}


Foam::tmp<Foam::pointField>
Foam::searchableSurfaceCollection::coordinates() const
{
    tmp<pointField> tCtrs = tmp<pointField>(new pointField(size()));
    pointField& ctrs = tCtrs();

    // Append individual coordinates
    label coordI = 0;

    forAll(subGeom_, surfI)
    {
        const pointField subCoords(subGeom_[surfI].coordinates());

        forAll(subCoords, i)
        {
            ctrs[coordI++] = transform_[surfI].globalPosition
            (
                cmptMultiply
                (
                    subCoords[i],
                    scale_[surfI]
                )
            );
        }
    }

    return tCtrs;
}


void Foam::searchableSurfaceCollection::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres.setSize(size());
    radiusSqr.setSize(centres.size());

    // Append individual coordinates
    label coordI = 0;

    forAll(subGeom_, surfI)
    {
        scalar maxScale = cmptMax(scale_[surfI]);

        pointField subCentres;
        scalarField subRadiusSqr;
        subGeom_[surfI].boundingSpheres(subCentres, subRadiusSqr);

        forAll(subCentres, i)
        {
            centres[coordI] = transform_[surfI].globalPosition
            (
                cmptMultiply
                (
                    subCentres[i],
                    scale_[surfI]
                )
            );
            radiusSqr[coordI] = maxScale*subRadiusSqr[i];
            coordI++;
        }
    }
}


Foam::tmp<Foam::pointField>
Foam::searchableSurfaceCollection::points() const
{
    // Get overall size
    label nPoints = 0;

    forAll(subGeom_, surfI)
    {
        nPoints += subGeom_[surfI].points()().size();
    }

    tmp<pointField> tPts(new pointField(nPoints));
    pointField& pts = tPts();

    // Append individual coordinates
    nPoints = 0;

    forAll(subGeom_, surfI)
    {
        const pointField subCoords(subGeom_[surfI].points());

        forAll(subCoords, i)
        {
            pts[nPoints++] = transform_[surfI].globalPosition
            (
                cmptMultiply
                (
                    subCoords[i],
                    scale_[surfI]
                )
            );
        }
    }

    return tPts;
}


void Foam::searchableSurfaceCollection::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& nearestInfo
) const
{
    // How to scale distance?
    scalarField minDistSqr(nearestDistSqr);

    labelList nearestSurf;
    findNearest
    (
        samples,
        minDistSqr,
        nearestInfo,
        nearestSurf
    );
}


void Foam::searchableSurfaceCollection::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());
    info = pointIndexHit();

    // Current nearest (to start) intersection
    pointField nearest(end);

    List<pointIndexHit> hitInfo(start.size());

    forAll(subGeom_, surfI)
    {
        // Starting point
        tmp<pointField> e0 = cmptDivide
        (
            transform_[surfI].localPosition
            (
                start
            ),
            scale_[surfI]
        );

        // Current best end point
        tmp<pointField> e1 = cmptDivide
        (
            transform_[surfI].localPosition
            (
                nearest
            ),
            scale_[surfI]
        );

        subGeom_[surfI].findLine(e0, e1, hitInfo);

        forAll(hitInfo, pointI)
        {
            if (hitInfo[pointI].hit())
            {
                // Transform back to global coordinate sys.
                nearest[pointI] = transform_[surfI].globalPosition
                (
                    cmptMultiply
                    (
                        hitInfo[pointI].rawPoint(),
                        scale_[surfI]
                    )
                );
                info[pointI] = hitInfo[pointI];
                info[pointI].rawPoint() = nearest[pointI];
                info[pointI].setIndex
                (
                    hitInfo[pointI].index()
                  + indexOffset_[surfI]
                );
            }
        }
    }


    // Debug check
    if (false)
    {
        forAll(info, pointI)
        {
            if (info[pointI].hit())
            {
                vector n(end[pointI] - start[pointI]);
                scalar magN = mag(n);

                if (magN > SMALL)
                {
                    n /= mag(n);

                    scalar s = ((info[pointI].rawPoint()-start[pointI])&n);

                    if (s < 0 || s > 1)
                    {
                        FatalErrorIn
                        (
                            "searchableSurfaceCollection::findLine(..)"
                        )   << "point:" << info[pointI]
                            << " s:" << s
                            << " outside vector "
                            << " start:" << start[pointI]
                            << " end:" << end[pointI]
                            << abort(FatalError);
                    }
                }
            }
        }
    }
}


void Foam::searchableSurfaceCollection::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    // To be done ...
    findLine(start, end, info);
}


void Foam::searchableSurfaceCollection::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit> >& info
) const
{
    // To be done. Assume for now only one intersection.
    List<pointIndexHit> nearestInfo;
    findLine(start, end, nearestInfo);

    info.setSize(start.size());
    forAll(info, pointI)
    {
        if (nearestInfo[pointI].hit())
        {
            info[pointI].setSize(1);
            info[pointI][0] = nearestInfo[pointI];
        }
        else
        {
            info[pointI].clear();
        }
    }
}


void Foam::searchableSurfaceCollection::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    if (subGeom_.size() == 0)
    {}
    else if (subGeom_.size() == 1)
    {
        if (mergeSubRegions_)
        {
            region.setSize(info.size());
            region = regionOffset_[0];
        }
        else
        {
            subGeom_[0].getRegion(info, region);
        }
    }
    else
    {
        // Multiple surfaces. Sort by surface.

        // Per surface the hit
        List<List<pointIndexHit> > surfInfo;
        // Per surface the original position
        List<List<label> > infoMap;
        sortHits(info, surfInfo, infoMap);

        region.setSize(info.size());
        region = -1;

        // Do region tests

        if (mergeSubRegions_)
        {
            // Actually no need for surfInfo. Just take region for surface.
            forAll(infoMap, surfI)
            {
                const labelList& map = infoMap[surfI];
                forAll(map, i)
                {
                    region[map[i]] = regionOffset_[surfI];
                }
            }
        }
        else
        {
            forAll(infoMap, surfI)
            {
                labelList surfRegion;
                subGeom_[surfI].getRegion(surfInfo[surfI], surfRegion);

                const labelList& map = infoMap[surfI];
                forAll(map, i)
                {
                    region[map[i]] = regionOffset_[surfI] + surfRegion[i];
                }
            }
        }
    }
}


void Foam::searchableSurfaceCollection::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    if (subGeom_.size() == 0)
    {}
    else if (subGeom_.size() == 1)
    {
        subGeom_[0].getNormal(info, normal);
    }
    else
    {
        // Multiple surfaces. Sort by surface.

        // Per surface the hit
        List<List<pointIndexHit> > surfInfo;
        // Per surface the original position
        List<List<label> > infoMap;
        sortHits(info, surfInfo, infoMap);

        normal.setSize(info.size());

        // Do region tests
        forAll(surfInfo, surfI)
        {
            vectorField surfNormal;
            subGeom_[surfI].getNormal(surfInfo[surfI], surfNormal);

            // Transform back to global coordinate sys.
            surfNormal = transform_[surfI].globalVector(surfNormal);

            const labelList& map = infoMap[surfI];
            forAll(map, i)
            {
                normal[map[i]] = surfNormal[i];
            }
        }
    }
}


void Foam::searchableSurfaceCollection::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    FatalErrorIn
    (
        "searchableSurfaceCollection::getVolumeType(const pointField&"
        ", List<volumeType>&) const"
    )   << "Volume type not supported for collection."
        << exit(FatalError);
}


void Foam::searchableSurfaceCollection::distribute
(
    const List<treeBoundBox>& bbs,
    const bool keepNonLocal,
    autoPtr<mapDistribute>& faceMap,
    autoPtr<mapDistribute>& pointMap
)
{
    forAll(subGeom_, surfI)
    {
        // Note:Tranform the bounding boxes? Something like
        // pointField bbPoints =
        // cmptDivide
        // (
        //     transform_[surfI].localPosition
        //     (
        //         bbs[i].points()
        //     ),
        //     scale_[surfI]
        // );
        // treeBoundBox newBb(bbPoints);

        // Note: what to do with faceMap, pointMap from multiple surfaces?
        subGeom_[surfI].distribute
        (
            bbs,
            keepNonLocal,
            faceMap,
            pointMap
        );
    }
}


void Foam::searchableSurfaceCollection::setField(const labelList& values)
{
    forAll(subGeom_, surfI)
    {
        subGeom_[surfI].setField
        (
            static_cast<const labelList&>
            (
                SubList<label>
                (
                    values,
                    subGeom_[surfI].size(),
                    indexOffset_[surfI]
                )
            )
        );
    }
}


void Foam::searchableSurfaceCollection::getField
(
    const List<pointIndexHit>& info,
    labelList& values
) const
{
    if (subGeom_.size() == 0)
    {}
    else if (subGeom_.size() == 1)
    {
        subGeom_[0].getField(info, values);
    }
    else
    {
        // Multiple surfaces. Sort by surface.

        // Per surface the hit
        List<List<pointIndexHit> > surfInfo;
        // Per surface the original position
        List<List<label> > infoMap;
        sortHits(info, surfInfo, infoMap);

        // Do surface tests
        forAll(surfInfo, surfI)
        {
            labelList surfValues;
            subGeom_[surfI].getField(surfInfo[surfI], surfValues);

            if (surfValues.size())
            {
                // Size values only when we have a surface that supports it.
                values.setSize(info.size());

                const labelList& map = infoMap[surfI];
                forAll(map, i)
                {
                    values[map[i]] = surfValues[i];
                }
            }
        }
    }
}


// ************************************************************************* //
