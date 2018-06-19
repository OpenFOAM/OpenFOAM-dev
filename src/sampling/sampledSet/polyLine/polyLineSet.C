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

#include "polyLineSet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyLineSet, 0);
    addToRunTimeSelectionTable(sampledSet, polyLineSet, word);

    const scalar polyLineSet::tol = 1e-6;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::polyLineSet::trackToBoundary
(
    passiveParticleCloud& particleCloud,
    passiveParticle& singleParticle,
    label& sampleI,
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<scalar>& samplingCurveDist
) const
{
    while (true)
    {
        // Local geometry info
        const vector offset = sampleCoords_[sampleI+1] - sampleCoords_[sampleI];
        const scalar smallDist = mag(tol*offset);

        singleParticle.track(offset, 0);
        const point trackPt = singleParticle.position();

        if (singleParticle.onBoundaryFace())
        {
            // Info<< "trackToBoundary : reached boundary"
            //    << "  trackPt:" << trackPt << endl;
            if
            (
                mag(trackPt - sampleCoords_[sampleI+1])
              < smallDist
            )
            {
                // Reached samplePt on boundary
                // Info<< "trackToBoundary : boundary. also sampling."
                //    << "  trackPt:" << trackPt << " sampleI+1:" << sampleI+1
                //    << endl;
                samplingPts.append(trackPt);
                samplingCells.append(singleParticle.cell());
                samplingFaces.append(singleParticle.face());

                // trackPt is at sampleI+1
                samplingCurveDist.append(1.0*(sampleI+1));
            }
            return true;
        }

        // Reached samplePt in cell
        samplingPts.append(trackPt);
        samplingCells.append(singleParticle.cell());
        samplingFaces.append(-1);

        // Convert trackPt to fraction in between sampleI and sampleI+1
        scalar dist =
            mag(trackPt - sampleCoords_[sampleI])
          / mag(sampleCoords_[sampleI+1] - sampleCoords_[sampleI]);
        samplingCurveDist.append(sampleI + dist);

        // go to next samplePt
        sampleI++;

        if (sampleI == sampleCoords_.size() - 1)
        {
            // no more samples.
            // Info<< "trackToBoundary : Reached end : sampleI now:" << sampleI
            //    << endl;
            return false;
        }
    }
}


void Foam::polyLineSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    // Check sampling points
    if (sampleCoords_.size() < 2)
    {
        FatalErrorInFunction
            << "Incorrect sample specification. Too few points:"
            << sampleCoords_ << exit(FatalError);
    }
    point oldPoint = sampleCoords_[0];
    for (label sampleI = 1; sampleI < sampleCoords_.size(); sampleI++)
    {
        if (mag(sampleCoords_[sampleI] - oldPoint) < small)
        {
            FatalErrorInFunction
                << "Incorrect sample specification."
                << " Point " << sampleCoords_[sampleI-1]
                << " at position " << sampleI-1
                << " and point " << sampleCoords_[sampleI]
                << " at position " << sampleI
                << " are too close" << exit(FatalError);
        }
        oldPoint = sampleCoords_[sampleI];
    }

    // Force calculation of cloud addressing on all processors
    const bool oldMoving = const_cast<polyMesh&>(mesh()).moving(false);
    passiveParticleCloud particleCloud(mesh());

    // current segment number
    label segmentI = 0;

    // starting index of current segment in samplePts
    label startSegmentI = 0;

    label sampleI = 0;

    point lastSample(great, great, great);
    while (true)
    {
        // Get boundary intersection
        point trackPt;
        label trackCelli = -1;
        label trackFacei = -1;

        do
        {
            const vector offset =
                sampleCoords_[sampleI+1] - sampleCoords_[sampleI];
            const scalar smallDist = mag(tol*offset);


            // Get all boundary intersections
            List<pointIndexHit> bHits = searchEngine().intersections
            (
                sampleCoords_[sampleI],
                sampleCoords_[sampleI+1]
            );

            point bPoint(great, great, great);
            label bFacei = -1;

            if (bHits.size())
            {
                bPoint = bHits[0].hitPoint();
                bFacei = bHits[0].index();
            }

            // Get tracking point

            bool isSample =
                getTrackingPoint
                (
                    sampleCoords_[sampleI],
                    bPoint,
                    bFacei,
                    smallDist,

                    trackPt,
                    trackCelli,
                    trackFacei
                );

            if (isSample && (mag(lastSample - trackPt) > smallDist))
            {
                // Info<< "calcSamples : getTrackingPoint returned valid sample"
                //     << "  trackPt:" << trackPt
                //     << "  trackFacei:" << trackFacei
                //     << "  trackCelli:" << trackCelli
                //     << "  sampleI:" << sampleI
                //     << "  dist:" << dist
                //     << endl;

                samplingPts.append(trackPt);
                samplingCells.append(trackCelli);
                samplingFaces.append(trackFacei);

                // Convert sampling position to unique curve parameter. Get
                // fraction of distance between sampleI and sampleI+1.
                scalar dist =
                    mag(trackPt - sampleCoords_[sampleI])
                  / mag(sampleCoords_[sampleI+1] - sampleCoords_[sampleI]);
                samplingCurveDist.append(sampleI + dist);

                lastSample = trackPt;
            }

            if (trackCelli == -1)
            {
                // No intersection found. Go to next point
                sampleI++;
            }
        } while ((trackCelli == -1) && (sampleI < sampleCoords_.size() - 1));

        if (sampleI == sampleCoords_.size() - 1)
        {
            // Info<< "calcSamples : Reached end of samples: "
            //    << "  sampleI now:" << sampleI
            //    << endl;
            break;
        }

        //
        // Segment sampleI .. sampleI+1 intersected by domain
        //

        // Initialize tracking starting from sampleI
        passiveParticle singleParticle
        (
            mesh(),
            trackPt,
            trackCelli
        );

        bool bReached = trackToBoundary
        (
            particleCloud,
            singleParticle,
            sampleI,
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

        if (!bReached)
        {
            // Info<< "calcSamples : Reached end of samples: "
            //    << "  sampleI now:" << sampleI
            //    << endl;
            break;
        }
        lastSample = singleParticle.position();


        // Find next boundary.
        sampleI++;

        if (sampleI == sampleCoords_.size() - 1)
        {
            // Info<< "calcSamples : Reached end of samples: "
            //    << "  sampleI now:" << sampleI
            //    << endl;
            break;
        }

        segmentI++;

        startSegmentI = samplingPts.size();
    }

    const_cast<polyMesh&>(mesh()).moving(oldMoving);
}


void Foam::polyLineSet::genSamples()
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

Foam::polyLineSet::polyLineSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const List<point>& sampleCoords
)
:
    sampledSet(name, mesh, searchEngine, axis),
    sampleCoords_(sampleCoords)
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


Foam::polyLineSet::polyLineSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    sampleCoords_(dict.lookup("points"))
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyLineSet::~polyLineSet()
{}


// ************************************************************************* //
