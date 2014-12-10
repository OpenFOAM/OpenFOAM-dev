/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "faceOnlySet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceOnlySet, 0);
    addToRunTimeSelectionTable(sampledSet, faceOnlySet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::faceOnlySet::trackToBoundary
(
    passiveParticleCloud& particleCloud,
    passiveParticle& singleParticle,
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<scalar>& samplingCurveDist
) const
{
    // distance vector between sampling points
    const vector offset = end_ - start_;
    const vector smallVec = tol*offset;
    const scalar smallDist = mag(smallVec);

    particle::TrackingData<passiveParticleCloud> trackData(particleCloud);

    // Alias
    const point& trackPt = singleParticle.position();

    while(true)
    {
        point oldPoint = trackPt;

        singleParticle.trackToFace(end_, trackData);

        if (singleParticle.face() != -1 && mag(oldPoint - trackPt) > smallDist)
        {
            // Reached face. Sample.
            samplingPts.append(trackPt);
            samplingCells.append(singleParticle.cell());
            samplingFaces.append(singleParticle.face());
            samplingCurveDist.append(mag(trackPt - start_));
        }

        if (mag(trackPt - end_) < smallDist)
        {
            // end reached
            return false;
        }
        else if (singleParticle.onBoundary())
        {
            // Boundary reached.
            return true;
        }
    }
}


void Foam::faceOnlySet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    // distance vector between sampling points
    if (mag(end_ - start_) < SMALL)
    {
        FatalErrorIn("faceOnlySet::calcSamples()")
            << "Incorrect sample specification :"
            << " start equals end point." << endl
            << "  start:" << start_
            << "  end:" << end_
            << exit(FatalError);
    }

    const vector offset = (end_ - start_);
    const vector normOffset = offset/mag(offset);
    const vector smallVec = tol*offset;
    const scalar smallDist = mag(smallVec);

    // Force calculation of cloud addressing on all processors
    const bool oldMoving = const_cast<polyMesh&>(mesh()).moving(false);
    passiveParticleCloud particleCloud(mesh());

    // Get all boundary intersections
    List<pointIndexHit> bHits = searchEngine().intersections
    (
        start_ - smallVec,
        end_ + smallVec
    );

    point bPoint(GREAT, GREAT, GREAT);
    label bFaceI = -1;

    if (bHits.size())
    {
        bPoint = bHits[0].hitPoint();
        bFaceI = bHits[0].index();
    }

    // Get first tracking point. Use bPoint, bFaceI if provided.

    point trackPt;
    label trackCellI = -1;
    label trackFaceI = -1;

    //Info<< "before getTrackingPoint : bPoint:" << bPoint
    //    << " bFaceI:" << bFaceI << endl;

    getTrackingPoint
    (
        offset,
        start_,
        bPoint,
        bFaceI,

        trackPt,
        trackCellI,
        trackFaceI
    );

    //Info<< "after getTrackingPoint : "
    //    << " trackPt:" << trackPt
    //    << " trackCellI:" << trackCellI
    //    << " trackFaceI:" << trackFaceI
    //    << endl;

    if (trackCellI == -1)
    {
        // Line start_ - end_ does not intersect domain at all.
        // (or is along edge)
        // Set points and cell/face labels to empty lists
        //Info<< "calcSamples : Both start_ and end_ outside domain"
        //    << endl;

        return;
    }

    if (trackFaceI == -1)
    {
        // No boundary face. Check for nearish internal face
        trackFaceI = findNearFace(trackCellI, trackPt, smallDist);
    }

    //Info<< "calcSamples : got first point to track from :"
    //    << "  trackPt:" << trackPt
    //    << "  trackCell:" << trackCellI
    //    << "  trackFace:" << trackFaceI
    //    << endl;

    //
    // Track until hit end of all boundary intersections
    //

    // current segment number
    label segmentI = 0;

    // starting index of current segment in samplePts
    label startSegmentI = 0;

    // index in bHits; current boundary intersection
    label bHitI = 1;

    while(true)
    {
        if (trackFaceI != -1)
        {
            //Info<< "trackPt:" << trackPt << " on face so use." << endl;
            samplingPts.append(trackPt);
            samplingCells.append(trackCellI);
            samplingFaces.append(trackFaceI);
            samplingCurveDist.append(mag(trackPt - start_));
        }

        // Initialize tracking starting from trackPt
        passiveParticle singleParticle
        (
            mesh(),
            trackPt,
            trackCellI
        );

        bool reachedBoundary = trackToBoundary
        (
            particleCloud,
            singleParticle,
            samplingPts,
            samplingCells,
            samplingFaces,
            samplingCurveDist
        );

        // fill sampleSegments
        for (label i = samplingPts.size() - 1; i >= startSegmentI; --i)
        {
            samplingSegments.append(segmentI);
        }


        if (!reachedBoundary)
        {
            //Info<< "calcSamples : Reached end of samples: "
            //    << "  samplePt now:" << singleParticle.position()
            //    << endl;
            break;
        }


        // Go past boundary intersection where tracking stopped
        // Use coordinate comparison instead of face comparison for
        // accuracy reasons

        bool foundValidB = false;

        while (bHitI < bHits.size())
        {
            scalar dist =
                (bHits[bHitI].hitPoint() - singleParticle.position())
              & normOffset;

            //Info<< "Finding next boundary : "
            //    << "bPoint:" << bHits[bHitI].hitPoint()
            //    << "  tracking:" << singleParticle.position()
            //    << "  dist:" << dist
            //    << endl;

            if (dist > smallDist)
            {
                // hitpoint is past tracking position
                foundValidB = true;
                break;
            }
            else
            {
                bHitI++;
            }
        }

        if (!foundValidB)
        {
            // No valid boundary intersection found beyond tracking position
            break;
        }

        // Update starting point for tracking
        trackFaceI = bHits[bHitI].index();
        trackPt = pushIn(bHits[bHitI].hitPoint(), trackFaceI);
        trackCellI = getBoundaryCell(trackFaceI);

        segmentI++;

        startSegmentI = samplingPts.size();
    }

    const_cast<polyMesh&>(mesh()).moving(oldMoving);
}


void Foam::faceOnlySet::genSamples()
{
    // Storage for sample points
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

    // Copy into *this
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

Foam::faceOnlySet::faceOnlySet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const point& start,
    const point& end
)
:
    sampledSet(name, mesh, searchEngine, axis),
    start_(start),
    end_(end)
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


Foam::faceOnlySet::faceOnlySet
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

Foam::faceOnlySet::~faceOnlySet()
{}


// ************************************************************************* //
