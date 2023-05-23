/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "refinementRegions.H"
#include "refinementSurfaces.H"
#include "searchableSurfaces.H"
#include "orientedSurface.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

template<>
const char*
NamedEnum<refinementRegions::refineMode, 5>::
names[] =
{
    "inside",
    "outside",
    "distance",
    "insideSpan",
    "outsideSpan"
};

}

const Foam::NamedEnum<Foam::refinementRegions::refineMode, 5>
    Foam::refinementRegions::refineModeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::refinementRegions::setAndCheckLevels
(
    const label shelli,
    const dictionary& dict
)
{
    const searchableSurface& shell = allGeometry_[shells_[shelli]];

    if
    (
        modes_[shelli] == refineMode::inside
     || modes_[shelli] == refineMode::outside
    )
    {
        if (!allGeometry_[shells_[shelli]].hasVolumeType())
        {
            WarningInFunction
                << "Shell " << shell.name()
                << " is not closed so testing for '"
                << refineModeNames_[modes_[shelli]]
                << "' may fail." << endl;
        }

        distances_[shelli].setSize(0);
        levels_[shelli].setSize(1);

        if (dict.found("levels") && !(dict.found("level")))
        {
            // Support 'levels' for backward compatibility
            const List<Tuple2<scalar, label>> distLevels(dict.lookup("levels"));

            if (distLevels.size() != 1)
            {
                FatalErrorInFunction
                    << "For refinement mode "
                    << refineModeNames_[modes_[shelli]]
                    << " specify only one distance+level."
                    << " (its distance gets discarded)"
                    << exit(FatalError);
            }

            levels_[shelli][0] = distLevels[0].second();
        }
        else
        {
            if (dict.found("levels"))
            {
                IOWarningInFunction(dict)
                    << "Found both 'level' and 'levels' entries, using 'level'."
                    << endl;
            }

            levels_[shelli][0] = readLabel(dict.lookup("level"));
        }

        if (modes_[shelli] == refineMode::inside)
        {
            Info<< "Refinement level " << levels_[shelli][0]
                << " for all cells inside " << shell.name() << endl;
        }
        else
        {
            Info<< "Refinement level " << levels_[shelli][0]
                << " for all cells outside " << shell.name() << endl;
        }
    }
    else if
    (
        modes_[shelli] == refineMode::insideSpan
     || modes_[shelli] == refineMode::outsideSpan
    )
    {
        if (!allGeometry_[shells_[shelli]].hasVolumeType())
        {
            WarningInFunction
                << "Shell " << shell.name()
                << " is not closed so testing for '"
                << refineModeNames_[modes_[shelli]]
                << "' may fail." << endl;
        }

        distances_[shelli].setSize(1);
        levels_[shelli].setSize(1);
        const Tuple2<scalar, label> distLevel(dict.lookup("level"));
        distances_[shelli][0] = distLevel.first();
        levels_[shelli][0] = distLevel.second();

        if (modes_[shelli] == refineMode::insideSpan)
        {
            Info<< "Refinement level " << levels_[shelli][0]
                << " for all cells inside " << shell.name()
                << " within distance " << distances_[shelli][0] << endl;
        }
        else
        {
            Info<< "Refinement level " << levels_[shelli][0]
                << " for all cells outside " << shell.name()
                << " within distance " << distances_[shelli][0] << endl;
        }
    }
    else
    {
        const List<Tuple2<scalar, label>> distLevels(dict.lookup("levels"));

        // Extract information into separate distance and level
        distances_[shelli].setSize(distLevels.size());
        levels_[shelli].setSize(distLevels.size());

        forAll(distLevels, j)
        {
            distances_[shelli][j] = distLevels[j].first();
            levels_[shelli][j] = distLevels[j].second();

            // Check in incremental order
            if (j > 0)
            {
                if
                (
                    (distances_[shelli][j] <= distances_[shelli][j-1])
                 || (levels_[shelli][j] > levels_[shelli][j-1])
                )
                {
                    FatalErrorInFunction
                        << "For refinement mode "
                        << refineModeNames_[modes_[shelli]]
                        << " : Refinement should be specified in order"
                        << " of increasing distance"
                        << " (and decreasing refinement level)." << endl
                        << "Distance:" << distances_[shelli][j]
                        << " refinementLevel:" << levels_[shelli][j]
                        << exit(FatalError);
                }
            }
        }

        if (modes_[shelli] == refineMode::distance)
        {
            Info<< "Refinement level according to distance to "
                << shell.name() << endl;

            forAll(levels_[shelli], j)
            {
                Info<< "    level " << levels_[shelli][j]
                    << " for all cells within " << distances_[shelli][j]
                    << " metre." << endl;
            }
        }
    }
}


void Foam::refinementRegions::orient()
{
    // Determine outside point.
    boundBox overallBb = boundBox::invertedBox;

    bool hasSurface = false;

    forAll(shells_, shelli)
    {
        const searchableSurface& s = allGeometry_[shells_[shelli]];

        if (modes_[shelli] != refineMode::distance && isA<triSurfaceMesh>(s))
        {
            const triSurfaceMesh& shell = refCast<const triSurfaceMesh>(s);

            if (shell.triSurface::size())
            {
                tmp<pointField> tpoints(shell.points());
                const pointField& points = tpoints();

                hasSurface = true;
                boundBox shellBb(points[0], points[0]);

                // Assume surface is compact!
                forAll(points, i)
                {
                    const point& pt = points[i];
                    shellBb.min() = min(shellBb.min(), pt);
                    shellBb.max() = max(shellBb.max(), pt);
                }

                overallBb.min() = min(overallBb.min(), shellBb.min());
                overallBb.max() = max(overallBb.max(), shellBb.max());
            }
        }
    }

    if (hasSurface)
    {
        const point outsidePt = overallBb.max() + overallBb.span();

        // Info<< "Using point " << outsidePt << " to orient shells" << endl;

        forAll(shells_, shelli)
        {
            const searchableSurface& s = allGeometry_[shells_[shelli]];

            if
            (
                modes_[shelli] != refineMode::distance
             && isA<triSurfaceMesh>(s)
            )
            {
                triSurfaceMesh& shell = const_cast<triSurfaceMesh&>
                (
                    refCast<const triSurfaceMesh>(s)
                );

                // Flip surface so outsidePt is outside.
                bool anyFlipped = orientedSurface::orient
                (
                    shell,
                    outsidePt,
                    true
                );

                if (anyFlipped)
                {
                    // orientedSurface will have done a clearOut of the surface.
                    // we could do a clearout of the triSurfaceMeshes::trees()
                    // but these aren't affected by orientation
                    // (except for cached
                    // sideness which should not be set at this point.
                    // !!Should check!)

                    Info<< "refinementRegions : Flipped orientation of surface "
                        << s.name()
                        << " so point " << outsidePt << " is outside." << endl;
                }
            }
        }
    }
}


Foam::scalar Foam::refinementRegions::interpolate
(
    const triSurfaceMesh& tsm,
    const triSurfacePointScalarField& closeness,
    const point& pt,
    const label index
) const
{
    const barycentric2D bary
    (
        triPointRef
        (
            tsm.points(),
            tsm.triSurface::operator[](index)
        ).pointToBarycentric(pt)
    );

    const labelledTri& lf = tsm.localFaces()[index];

    return closeness[lf[0]]*bary[0]
         + closeness[lf[1]]*bary[1]
         + closeness[lf[2]]*bary[2];
}


void Foam::refinementRegions::findHigherLevel
(
    const pointField& pt,
    const label shelli,
    const scalar level0EdgeLength,
    labelList& maxLevel
) const
{
    const labelList& levels = levels_[shelli];

    if (modes_[shelli] == refineMode::distance)
    {
        // Distance mode.

        const scalarField& distances = distances_[shelli];

        // Collect all those points that have a current maxLevel less than
        // (any of) the shell. Also collect the furthest distance allowable
        // to any shell with a higher level.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        scalarField candidateDistSqr(pt.size());
        label candidatei = 0;

        forAll(maxLevel, pointi)
        {
            forAllReverse(levels, leveli)
            {
                if (levels[leveli] > maxLevel[pointi])
                {
                    candidates[candidatei] = pt[pointi];
                    candidateMap[candidatei] = pointi;
                    candidateDistSqr[candidatei] = sqr(distances[leveli]);
                    candidatei++;
                    break;
                }
            }
        }
        candidates.setSize(candidatei);
        candidateMap.setSize(candidatei);
        candidateDistSqr.setSize(candidatei);

        // Do the expensive nearest test only for the candidate points.
        List<pointIndexHit> nearInfo;
        allGeometry_[shells_[shelli]].findNearest
        (
            candidates,
            candidateDistSqr,
            nearInfo
        );

        // Update maxLevel
        forAll(nearInfo, candidatei)
        {
            if (nearInfo[candidatei].hit())
            {
                // Check which level it actually is in.
                label minDistI = findLower
                (
                    distances,
                    mag(nearInfo[candidatei].hitPoint()-candidates[candidatei])
                );

                label pointi = candidateMap[candidatei];

                // pt is in between shell[minDistI] and shell[minDistI+1]
                maxLevel[pointi] = levels[minDistI+1];
            }
        }
    }
    else if
    (
        modes_[shelli] == refineMode::insideSpan
     || modes_[shelli] == refineMode::outsideSpan
    )
    {
        const triSurfaceMesh& tsm =
            refCast<const triSurfaceMesh>(allGeometry_[shells_[shelli]]);

        // Collect all those points that have a current maxLevel less than
        // the maximum and furthest distance allowable for the shell.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        scalarField candidateDistSqr(pt.size());
        label candidatei = 0;

        forAll(pt, pointi)
        {
            if (levels[0] > maxLevel[pointi])
            {
                candidates[candidatei] = pt[pointi];
                candidateMap[candidatei] = pointi;
                candidateDistSqr[candidatei] = sqr(distances_[shelli][0]);
                candidatei++;
            }
        }
        candidates.setSize(candidatei);
        candidateMap.setSize(candidatei);
        candidateDistSqr.setSize(candidatei);

        // Do the expensive nearest test only for the candidate points.
        List<pointIndexHit> nearInfo;
        tsm.findNearest
        (
            candidates,
            candidateDistSqr,
            nearInfo
        );

        // Minimum span for the maximum level specified
        const scalar minSpan
        (
            cellsAcrossSpan_[shelli]*level0EdgeLength/(1 << (levels[0] - 1))
        );

        // Update maxLevel
        forAll(nearInfo, candidatei)
        {
            if (nearInfo[candidatei].hit())
            {
                const scalar span
                (
                    interpolate
                    (
                        tsm,
                        closeness_[shelli],
                        nearInfo[candidatei].rawPoint(),
                        nearInfo[candidatei].index()
                    )
                );

                if (span > minSpan)
                {
                    const label level
                    (
                        log2(cellsAcrossSpan_[shelli]*level0EdgeLength/span) + 1
                    );

                    maxLevel[candidateMap[candidatei]] = min(levels[0], level);
                }
                else
                {
                    maxLevel[candidateMap[candidatei]] = levels[0];
                }
            }
        }
    }
    else
    {
        // Inside/outside mode

        // Collect all those points that have a current maxLevel less than the
        // shell.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        label candidatei = 0;

        forAll(maxLevel, pointi)
        {
            if (levels[0] > maxLevel[pointi])
            {
                candidates[candidatei] = pt[pointi];
                candidateMap[candidatei] = pointi;
                candidatei++;
            }
        }
        candidates.setSize(candidatei);
        candidateMap.setSize(candidatei);

        // Do the expensive nearest test only for the candidate points.
        List<volumeType> volType;
        allGeometry_[shells_[shelli]].getVolumeType(candidates, volType);

        forAll(volType, i)
        {
            label pointi = candidateMap[i];

            if
            (
                (
                    modes_[shelli] == refineMode::inside
                 && volType[i] == volumeType::inside
                )
             || (
                    modes_[shelli] == refineMode::outside
                 && volType[i] == volumeType::outside
                )
            )
            {
                maxLevel[pointi] = levels[0];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementRegions::refinementRegions
(
    const searchableSurfaces& allGeometry,
    const dictionary& shellsDict
)
:
    allGeometry_(allGeometry)
{
    // Wildcard specification : loop over all surfaces and try to find a match.

    // Count number of shells.
    label shelli = 0;
    forAll(allGeometry.names(), geomi)
    {
        const word& geomName = allGeometry_.names()[geomi];

        if (shellsDict.found(geomName))
        {
            shelli++;
        }
    }

    // Size lists
    shells_.setSize(shelli);
    modes_.setSize(shelli);
    distances_.setSize(shelli);
    levels_.setSize(shelli);
    cellsAcrossSpan_.setSize(shelli);
    closeness_.setSize(shelli);

    HashSet<word> unmatchedKeys(shellsDict.toc());
    shelli = 0;

    forAll(allGeometry_.names(), geomi)
    {
        const word& geomName = allGeometry_.names()[geomi];

        const entry* ePtr = shellsDict.lookupEntryPtr(geomName, false, true);

        if (ePtr)
        {
            const dictionary& dict = ePtr->dict();
            unmatchedKeys.erase(ePtr->keyword());

            shells_[shelli] = geomi;
            modes_[shelli] = refineModeNames_.read(dict.lookup("mode"));

            // Read pairs of distance+level
            setAndCheckLevels(shelli, dict);

            if
            (
                modes_[shelli] == refineMode::insideSpan
             || modes_[shelli] == refineMode::outsideSpan
            )
            {
                const searchableSurface& surface = allGeometry_[geomi];

                if (isA<triSurfaceMesh>(surface))
                {
                    dict.lookup("cellsAcrossSpan") >> cellsAcrossSpan_[shelli];

                    const triSurfaceMesh& tsm =
                        refCast<const triSurfaceMesh>(surface);

                    closeness_.set
                    (
                        shelli,
                        new triSurfacePointScalarField
                        (
                            IOobject
                            (
                                surface.name()
                              + ".closeness."
                              + (
                                    modes_[shelli] == refineMode::insideSpan
                                  ? "internal"
                                  : "external"
                                )
                              + "PointCloseness",
                                surface.searchableSurface::time().constant(),
                                searchableSurface::geometryDir
                                (
                                    surface.searchableSurface::time()
                                ),
                                surface.searchableSurface::time(),
                                IOobject::MUST_READ
                            ),
                            tsm,
                            dimLength,
                            true
                        )
                    );
                }
                else
                {
                    FatalIOErrorInFunction(shellsDict)
                        << "Surface " << surface.name()
                        << " is not a triSurface as required by"
                           " refinement modes "
                        << refineModeNames_[refineMode::insideSpan]
                        << " and " << refineModeNames_[refineMode::outsideSpan]
                        << exit(FatalIOError);
                }
            }

            shelli++;
        }
    }

    if (unmatchedKeys.size() > 0)
    {
        IOWarningInFunction
        (
            shellsDict
        )   << "Not all entries in refinementRegions dictionary were used."
            << " The following entries were not used : "
            << unmatchedKeys.sortedToc()
            << endl;
    }

    // Orient shell surfaces before any searching is done. Note that this
    // only needs to be done for inside or outside. Orienting surfaces
    // constructs lots of addressing which we want to avoid.
    orient();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::refinementRegions::maxLevel() const
{
    label overallMax = 0;

    forAll(levels_, shelli)
    {
        overallMax = max(overallMax, max(levels_[shelli]));
    }

    return overallMax;
}


void Foam::refinementRegions::findHigherLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    const scalar level0EdgeLength,
    labelList& maxLevel
) const
{
    // Maximum level of any shell. Start off with level of point.
    maxLevel = ptLevel;

    forAll(shells_, shelli)
    {
        findHigherLevel(pt, shelli, level0EdgeLength, maxLevel);
    }
}


// ************************************************************************* //
