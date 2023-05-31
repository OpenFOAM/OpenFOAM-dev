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

#include "lineFace.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "treeDataCell.H"
#include "sampledSetCloud.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{
    defineTypeNameAndDebug(lineFace, 0);
    addToRunTimeSelectionTable(sampledSet, lineFace, word);
}
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::sampledSets::lineFace::calcSamples
(
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const vector& start,
    const vector& end,
    const label storeFaces,
    const bool storeCells,
    DynamicList<point>& samplingPositions,
    DynamicList<scalar>& samplingDistances,
    DynamicList<label>& samplingSegments,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces
)
{
    // Get all candidates for starting the tracks
    List<DynamicList<label>> procCandidateCells(Pstream::nProcs());
    List<DynamicList<scalar>> procCandidateTs(Pstream::nProcs());
    {
        const label startCelli = searchEngine.findCell(start);
        if (startCelli != -1)
        {
            procCandidateCells[Pstream::myProcNo()].append(startCelli);
            procCandidateTs[Pstream::myProcNo()].append(0);
        }

        const label endCelli = searchEngine.findCell(end);
        if (endCelli != -1)
        {
            procCandidateCells[Pstream::myProcNo()].append(endCelli);
            procCandidateTs[Pstream::myProcNo()].append(1);
        }

        const List<pointIndexHit> bHits =
            searchEngine.intersections(start, end);
        forAll(bHits, bHiti)
        {
            for (label bHitj = bHiti + 1; bHitj < bHits.size(); ++ bHitj)
            {
                const point midP =
                    (bHits[bHiti].hitPoint() + bHits[bHitj].hitPoint())/2;

                const label midCelli = searchEngine.findCell(midP);
                const scalar midT = mag(midP - start)/mag(end - start);

                if (midCelli != -1)
                {
                    procCandidateCells[Pstream::myProcNo()].append(midCelli);
                    procCandidateTs[Pstream::myProcNo()].append(midT);
                }
            }
        }
    }
    Pstream::gatherList(procCandidateCells);
    Pstream::scatterList(procCandidateCells);
    Pstream::gatherList(procCandidateTs);
    Pstream::scatterList(procCandidateTs);

    // Tracking data
    const List<point> startEnd({start, end});
    const List<point> endStart({end, start});
    DynamicList<point> halfSegmentPositions;
    DynamicList<scalar> halfSegmentDistances;
    DynamicList<label> halfSegmentCells;
    DynamicList<label> halfSegmentFaces;

    // Create a cloud with which to track segments
    sampledSetCloud particles
    (
        mesh,
        lineFace::typeName,
        IDLList<sampledSetParticle>()
    );

    // Create each segment in turn
    label segmenti = 0;
    forAll(procCandidateCells, proci)
    {
        forAll(procCandidateCells[proci], candidatei)
        {
            const label celli = procCandidateCells[proci][candidatei];

            if (celli == -1) continue;

            const scalar t = procCandidateTs[proci][candidatei];
            const point p = (1 - t)*start + t*end;

            // Track in both directions to form parts of the segment either
            // side of the candidate
            Pair<scalar> segmentT;
            forAll(segmentT, i)
            {
                sampledSetParticle::trackingData tdBwd
                (
                    particles,
                    i == 0 ? endStart : startEnd,
                    false,
                    storeFaces,
                    storeCells,
                    halfSegmentPositions,
                    halfSegmentDistances,
                    halfSegmentCells,
                    halfSegmentFaces
                );

                particles.clear();
                if (proci == Pstream::myProcNo())
                {
                    particles.addParticle
                    (
                        new sampledSetParticle
                        (
                            mesh,
                            p,
                            celli,
                            0,
                            i == 0 ? t : 1 - t,
                            0
                        )
                    );
                }

                particles.move(particles, tdBwd);

                segmentT[i] =
                    returnReduce
                    (
                        particles.size()
                      ? mag(particles.first()->position(mesh) - start)
                       /mag(end - start)
                      : i == 0 ? vGreat : -vGreat,
                        [i](const scalar a, const scalar b)
                        {
                            return (i == 0) == (a < b) ? a : b;
                        }
                    );

                if (i == 0)
                {
                    reverse(halfSegmentPositions);
                    reverse(halfSegmentDistances);
                    reverse(halfSegmentCells);
                    reverse(halfSegmentFaces);
                }

                const label n = halfSegmentPositions.size();

                samplingPositions.append(halfSegmentPositions);
                samplingSegments.append(labelList(n, segmenti));
                samplingCells.append(halfSegmentCells);
                samplingFaces.append(halfSegmentFaces);

                halfSegmentPositions.clear();
                halfSegmentDistances.clear();
                halfSegmentCells.clear();
                halfSegmentFaces.clear();

                // If storing cells we need to store the starting cells between
                // the tracks
                if (proci == Pstream::myProcNo() && i == 0 && storeCells)
                {
                    particle trackBwd(mesh, p, celli), trackFwd(trackBwd);
                    trackBwd.trackToFace(mesh, start - p, 0);
                    trackFwd.trackToFace(mesh, end - p, 0);
                    if (trackBwd.onFace() && trackFwd.onFace())
                    {
                        samplingPositions.append
                        (
                            (
                                trackBwd.position(mesh)
                              + trackFwd.position(mesh)
                            )/2
                        );
                        samplingSegments.append(segmenti);
                        samplingCells.append(celli);
                        samplingFaces.append(-1);
                    }
                }
            }

            // Disable all candidates that fall within the bounds of the
            // computed segment
            forAll(procCandidateCells, procj)
            {
                forAll(procCandidateCells[procj], candidatej)
                {
                    const label cellj = procCandidateCells[procj][candidatej];

                    if (cellj == -1) continue;

                    const scalar t = procCandidateTs[procj][candidatej];

                    if
                    (
                        t > segmentT.first() - rootSmall
                     && t < segmentT.second() + rootSmall
                    )
                    {
                        procCandidateCells[procj][candidatej] = -1;
                        procCandidateTs[procj][candidatej] = NaN;
                    }
                }
            }

            // Move onto the next segment
            segmenti ++;
        }
    }

    // Set the distances
    samplingDistances.resize(samplingPositions.size());
    forAll(samplingPositions, i)
    {
        samplingDistances[i] = mag(samplingPositions[i] - start);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::lineFace::calcSamples
(
    DynamicList<point>& samplingPositions,
    DynamicList<scalar>& samplingDistances,
    DynamicList<label>& samplingSegments,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces
) const
{
    calcSamples
    (
        mesh(),
        searchEngine(),
        start_,
        end_,
        1,
        false,
        samplingPositions,
        samplingDistances,
        samplingSegments,
        samplingCells,
        samplingFaces
    );
}


void Foam::sampledSets::lineFace::genSamples()
{
    DynamicList<point> samplingPositions;
    DynamicList<scalar> samplingDistances;
    DynamicList<label> samplingSegments;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;

    calcSamples
    (
        samplingPositions,
        samplingDistances,
        samplingSegments,
        samplingCells,
        samplingFaces
    );

    samplingPositions.shrink();
    samplingDistances.shrink();
    samplingSegments.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();

    setSamples
    (
        samplingPositions,
        samplingDistances,
        samplingSegments,
        samplingCells,
        samplingFaces
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSets::lineFace::lineFace
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    start_(dict.lookup("start")),
    end_(dict.lookup("end"))
{
    genSamples();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::lineFace::~lineFace()
{}


// ************************************************************************* //
