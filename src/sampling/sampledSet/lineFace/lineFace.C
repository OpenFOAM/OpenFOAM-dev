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
    // Look for a starting location within the mesh
    vector startPt = start;
    label startFaceI = -1, startCellI = searchEngine.findCell(start);

    // Get all the hits of the line with the boundaries
    const List<pointIndexHit> bHits =
        searchEngine.intersections(start, end);

    // Loop over the hits, starting new segments each time
    label bHitI = startCellI == -1 ? 0 : -1;
    label sampleSegmentI = 0;
    for (; bHitI < bHits.size(); ++ bHitI)
    {
        // If a boundary start point, then initialise to the current hit
        if (bHitI >= 0)
        {
            startPt = bHits[bHitI].hitPoint();
            startFaceI = bHits[bHitI].index();
            startCellI = mesh.faceOwner()[startFaceI];
        }

        // If the hit points outward, move on to the next one
        if (startFaceI != -1)
        {
            const vector bNormal = normalised(mesh.faceAreas()[startFaceI]);
            const scalar bDot = bNormal & (end - start);
            if (bDot > 0)
            {
                continue;
            }
        }

        // Create a particle. If we are starting on a boundary face, track
        // backwards into it so that the particle has the correct topology.
        passiveParticle sampleParticle(mesh, startPt, startCellI);
        if (startFaceI != -1)
        {
            sampleParticle.track(start - end, 0);
            if (!sampleParticle.onBoundaryFace())
            {
                FatalErrorInFunction
                    << "Failed to associate with the starting boundary face"
                    << exit(FatalError);
            }
        }

        // Track until a boundary is hit, appending the face intersections to
        // the lists of samples
        while (true)
        {
            const point pt = sampleParticle.position();
            const scalar dist = mag(pt - start);
            const bool first =
                samplingSegments.size() == 0
             || samplingSegments.last() != sampleSegmentI;

            if (sampleParticle.onFace())
            {
                samplingPts.append(pt);
                samplingCells.append(sampleParticle.cell());
                samplingFaces.append(sampleParticle.face());
                samplingSegments.append(sampleSegmentI);
                samplingCurveDist.append(dist);
            }

            const vector s = (1 - dist/mag(end - start))*(end - start);

            if
            (
                (!first && sampleParticle.onBoundaryFace())
             || sampleParticle.trackToCell(s, 0) == 0
            )
            {
                break;
            }
        }

        ++ sampleSegmentI;
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
