/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "searchableSurfaceWithGaps.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(searchableSurfaceWithGaps, 0);
addToRunTimeSelectionTable(searchableSurface, searchableSurfaceWithGaps, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Pair<Foam::vector> Foam::searchableSurfaceWithGaps::offsetVecs
(
    const point& start,
    const point& end
) const
{
    Pair<vector> offsets(vector::zero, vector::zero);

    vector n(end-start);

    scalar magN = mag(n);

    if (magN > SMALL)
    {
        n /= magN;

        // Do first offset vector. Is the coordinate axes with the smallest
        // component along the vector n.
        scalar minMag = GREAT;
        direction minCmpt = 0;

        for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
        {
            if (mag(n[cmpt]) < minMag)
            {
                minMag = mag(n[cmpt]);
                minCmpt = cmpt;
            }
        }

        offsets[0][minCmpt] = 1.0;
        // Orthonormalise
        offsets[0] -= n[minCmpt]*n;
        offsets[0] /= mag(offsets[0]);
        // Do second offset vector perp to original edge and first offset vector
        offsets[1] = n ^ offsets[0];

        // Scale
        offsets[0] *= gap_;
        offsets[1] *= gap_;
    }

    return offsets;
}


void Foam::searchableSurfaceWithGaps::offsetVecs
(
    const pointField& start,
    const pointField& end,
    pointField& offset0,
    pointField& offset1
) const
{
    offset0.setSize(start.size());
    offset1.setSize(start.size());

    forAll(start, i)
    {
        const Pair<vector> offsets(offsetVecs(start[i], end[i]));
        offset0[i] = offsets[0];
        offset1[i] = offsets[1];
    }
}


Foam::label Foam::searchableSurfaceWithGaps::countMisses
(
    const List<pointIndexHit>& info,
    labelList& missMap
)
{
    label nMiss = 0;
    forAll(info, i)
    {
        if (!info[i].hit())
        {
            nMiss++;
        }
    }

    missMap.setSize(nMiss);
    nMiss = 0;

    forAll(info, i)
    {
        if (!info[i].hit())
        {
            missMap[nMiss++] = i;
        }
    }

    return nMiss;
}


// Anything not a hit in both counts as a hit
Foam::label Foam::searchableSurfaceWithGaps::countMisses
(
    const List<pointIndexHit>& plusInfo,
    const List<pointIndexHit>& minInfo,
    labelList& missMap
)
{
    label nMiss = 0;
    forAll(plusInfo, i)
    {
        if (!plusInfo[i].hit() || !minInfo[i].hit())
        {
            nMiss++;
        }
    }

    missMap.setSize(nMiss);
    nMiss = 0;

    forAll(plusInfo, i)
    {
        if (!plusInfo[i].hit() || !minInfo[i].hit())
        {
            missMap[nMiss++] = i;
        }
    }

    return nMiss;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceWithGaps::searchableSurfaceWithGaps
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    gap_(readScalar(dict.lookup("gap"))),
    subGeom_(1)
{
    const word subGeomName(dict.lookup("surface"));

    const searchableSurface& s =
        io.db().lookupObject<searchableSurface>(subGeomName);

    subGeom_.set(0, &const_cast<searchableSurface&>(s));

    bounds() = subGeom_[0].bounds();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableSurfaceWithGaps::~searchableSurfaceWithGaps()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::searchableSurfaceWithGaps::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{

    // Test with unperturbed vectors
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    surface().findLine(start, end, info);

    // Count number of misses. Determine map
    labelList compactMap;
    label nMiss = countMisses(info, compactMap);

    if (returnReduce(nMiss, sumOp<label>()) > 0)
    {
        //Pout<< "** retesting with offset0 " << nMiss << " misses out of "
        //    << start.size() << endl;

        // extract segments according to map
        pointField compactStart(start, compactMap);
        pointField compactEnd(end, compactMap);

        // Calculate offset vector
        pointField offset0, offset1;
        offsetVecs
        (
            compactStart,
            compactEnd,
            offset0,
            offset1
        );

        // Test with offset0 perturbed vectors
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // test in pairs: only if both perturbations hit something
        // do we accept the hit.

        const vectorField smallVec(1e-6*(compactEnd-compactStart));

        List<pointIndexHit> plusInfo;
        surface().findLine
        (
            compactStart+offset0-smallVec,
            compactEnd+offset0+smallVec,
            plusInfo
        );
        List<pointIndexHit> minInfo;
        surface().findLine
        (
            compactStart-offset0-smallVec,
            compactEnd-offset0+smallVec,
            minInfo
        );

        // Extract any hits
        forAll(plusInfo, i)
        {
            if (plusInfo[i].hit() && minInfo[i].hit())
            {
                info[compactMap[i]] = plusInfo[i];
                info[compactMap[i]].rawPoint() -= offset0[i];
            }
        }

        labelList plusMissMap;
        nMiss = countMisses(plusInfo, minInfo, plusMissMap);

        if (returnReduce(nMiss, sumOp<label>()) > 0)
        {
            //Pout<< "** retesting with offset1 " << nMiss << " misses out of "
            //    << start.size() << endl;

            // Test with offset1 perturbed vectors
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            // Extract (inplace possible because of order)
            forAll(plusMissMap, i)
            {
                label mapI = plusMissMap[i];
                compactStart[i] = compactStart[mapI];
                compactEnd[i] = compactEnd[mapI];
                compactMap[i] = compactMap[mapI];
                offset0[i] = offset0[mapI];
                offset1[i] = offset1[mapI];
            }
            compactStart.setSize(plusMissMap.size());
            compactEnd.setSize(plusMissMap.size());
            compactMap.setSize(plusMissMap.size());
            offset0.setSize(plusMissMap.size());
            offset1.setSize(plusMissMap.size());

            const vectorField smallVec(1e-6*(compactEnd-compactStart));

            surface().findLine
            (
                compactStart+offset1-smallVec,
                compactEnd+offset1+smallVec,
                plusInfo
            );
            surface().findLine
            (
                compactStart-offset1-smallVec,
                compactEnd-offset1+smallVec,
                minInfo
            );

            // Extract any hits
            forAll(plusInfo, i)
            {
                if (plusInfo[i].hit() && minInfo[i].hit())
                {
                    info[compactMap[i]] = plusInfo[i];
                    info[compactMap[i]].rawPoint() -= offset1[i];
                }
            }
        }
    }
}


void Foam::searchableSurfaceWithGaps::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    // To be done ...
    findLine(start, end, info);
}


void Foam::searchableSurfaceWithGaps::findLineAll
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


// ************************************************************************* //
