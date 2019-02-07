/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
)
{
    // Create lists of initial positions from which to track, the faces and
    // cells associated with those positions, and whether the track  propagates
    // forward (true) or backward (false) along the line from start to end
    DynamicList<point> initialPts;
    DynamicList<label> initialFaces, initialCells;
    DynamicList<bool> initialDirections;

    // Add boundary hits
    const List<pointIndexHit> bHits = searchEngine.intersections(start, end);
    forAll(bHits, bHiti)
    {
        initialPts.append(bHits[bHiti].hitPoint());
        const label facei = bHits[bHiti].index();
        initialFaces.append(facei);
        initialCells.append(mesh.faceOwner()[facei]);
        initialDirections.append((mesh.faceAreas()[facei] & (end - start)) < 0);
    }

    // Add the start and end points if they can be found within the mesh
    const label startCelli = searchEngine.findCell(start);
    if (startCelli != -1)
    {
        initialPts.append(start);
        initialFaces.append(-1);
        initialCells.append(startCelli);
        initialDirections.append(true);
    }
    const label endCelli = searchEngine.findCell(end);
    if (endCelli != -1)
    {
        initialPts.append(end);
        initialFaces.append(-1);
        initialCells.append(endCelli);
        initialDirections.append(false);
    }

    // Loop over the initial points, starting new segments each time
    label sampleSegmenti = 0;
    DynamicList<Pair<point>> lines;
    forAll(initialPts, initiali)
    {
        // Get the sign
        const scalar sign = initialDirections[initiali] ? +1 : -1;

        // Create a particle. Track backwards into the boundary face so that
        // the particle has the correct topology.
        passiveParticle sampleParticle
        (
            mesh,
            initialPts[initiali],
            initialCells[initiali]
        );
        if (initialFaces[initiali] != -1)
        {
            sampleParticle.track(sign*(start - end), 0);
            if (!sampleParticle.onBoundaryFace())
            {
                FatalErrorInFunction
                    << "Failed to associate with the starting boundary face"
                    << exit(FatalError);
            }
        }

        // Track until a boundary is hit, appending the face intersections
        // to the lists of samples, and storing the line
        DynamicList<point> segmentPts;
        DynamicList<label> segmentCells, segmentFaces;
        Pair<point> line(sampleParticle.position(), sampleParticle.position());
        while (true)
        {
            const point pt = sampleParticle.position();
            const scalar dist = mag(pt - (sign > 0 ? start : end));
            const bool first = segmentPts.size() == 0;

            if (sampleParticle.onFace())
            {
                segmentPts.append(pt);
                segmentCells.append(sampleParticle.cell());
                segmentFaces.append(sampleParticle.face());
            }

            const vector s =
                sign*(end - start)*(1 - dist/mag(end - start));

            if
            (
                (!first && sampleParticle.onBoundaryFace())
             || sampleParticle.trackToCell(s, 0) == 0
            )
            {
                break;
            }
        }
        line[1] = sampleParticle.position();

        // Reverse if going backwards
        if (sign < 0)
        {
            inplaceReverseList(segmentPts);
            inplaceReverseList(segmentCells);
            inplaceReverseList(segmentFaces);
            line = reverse(line);
        }

        // Mark point as not to be kept if they fall within the bounds of
        // previous lines
        boolList segmentKeep(segmentPts.size(), true);
        forAll(segmentPts, segmentPti)
        {
            forAll(lines, linei)
            {
                const Pair<point>& l = lines[linei];
                const vector dlHat = normalised(l[1] - l[0]);
                if (magSqr(dlHat) == 0)
                {
                    continue;
                }
                const scalar dot0 = (segmentPts[segmentPti] - l[0]) & dlHat;
                const scalar dot1 = (l[1] - segmentPts[segmentPti]) & dlHat;
                if (dot0 > 0 && dot1 > 0)
                {
                    segmentKeep[segmentPti] = false;
                    break;
                }
            }
        }

        // Store the line
        lines.append(line);

        // Add new segments to the lists, breaking the segment anywhere that
        // points are not kept
        bool newSampleSegment = false;
        forAll(segmentPts, segmentPti)
        {
            if (segmentKeep[segmentPti])
            {
                samplingPts.append(segmentPts[segmentPti]);
                samplingCells.append(segmentCells[segmentPti]);
                samplingFaces.append(segmentFaces[segmentPti]);
                samplingSegments.append(sampleSegmenti);
                samplingCurveDist.append(mag(segmentPts[segmentPti] - start));
                newSampleSegment = true;
            }
            else if (newSampleSegment)
            {
                ++ sampleSegmenti;
                newSampleSegment = false;
            }
        }
        if (newSampleSegment)
        {
            ++ sampleSegmenti;
            newSampleSegment = false;
        }
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::lineFace::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    calcSamples
    (
        mesh(),
        searchEngine(),
        start_,
        end_,
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );
}


void Foam::sampledSets::lineFace::genSamples()
{
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    calcSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    setSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
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

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::lineFace::~lineFace()
{}


// ************************************************************************* //
