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

#include "sampledSet.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "meshSearch.H"
#include "writer.H"
#include "particle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    const scalar sampledSet::tol = 1e-6;

    defineTypeNameAndDebug(sampledSet, 0);
    defineRunTimeSelectionTable(sampledSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::sampledSet::getBoundaryCell(const label faceI) const
{
    return mesh().faceOwner()[faceI];
}


Foam::label Foam::sampledSet::getCell
(
    const label faceI,
    const point& sample
) const
{
    if (faceI == -1)
    {
        FatalErrorIn
        (
            "sampledSet::getCell(const label, const point&)"
        )   << "Illegal face label " << faceI
            << abort(FatalError);
    }

    if (faceI >= mesh().nInternalFaces())
    {
        label cellI = getBoundaryCell(faceI);

        if (!mesh().pointInCell(sample, cellI, searchEngine_.decompMode()))
        {
            FatalErrorIn
            (
                "sampledSet::getCell(const label, const point&)"
            )   << "Found cell " << cellI << " using face " << faceI
                << ". But cell does not contain point " << sample
                << abort(FatalError);
        }
        return cellI;
    }
    else
    {
        // Try owner and neighbour to see which one contains sample

        label cellI = mesh().faceOwner()[faceI];

        if (mesh().pointInCell(sample, cellI, searchEngine_.decompMode()))
        {
            return cellI;
        }
        else
        {
            cellI = mesh().faceNeighbour()[faceI];

            if (mesh().pointInCell(sample, cellI, searchEngine_.decompMode()))
            {
                return cellI;
            }
            else
            {
                FatalErrorIn
                (
                    "sampledSet::getCell(const label, const point&)"
                )   << "None of the neighbours of face "
                    << faceI << " contains point " << sample
                    << abort(FatalError);

                return -1;
            }
        }
    }
}


Foam::scalar Foam::sampledSet::calcSign
(
    const label faceI,
    const point& sample
) const
{
    vector vec = sample - mesh().faceCentres()[faceI];

    scalar magVec = mag(vec);

    if (magVec < VSMALL)
    {
        // sample on face centre. Regard as inside
        return -1;
    }

    vec /= magVec;

    vector n = mesh().faceAreas()[faceI];

    n /= mag(n) + VSMALL;

    return n & vec;
}


// Return face (or -1) of face which is within smallDist of sample
Foam::label Foam::sampledSet::findNearFace
(
    const label cellI,
    const point& sample,
    const scalar smallDist
) const
{
    const cell& myFaces = mesh().cells()[cellI];

    forAll(myFaces, myFaceI)
    {
        const face& f = mesh().faces()[myFaces[myFaceI]];

        pointHit inter = f.nearestPoint(sample, mesh().points());

        scalar dist;

        if (inter.hit())
        {
            dist = mag(inter.hitPoint() - sample);
        }
        else
        {
            dist = mag(inter.missPoint() - sample);
        }

        if (dist < smallDist)
        {
            return myFaces[myFaceI];
        }
    }
    return -1;
}


// 'Pushes' point facePt (which is almost on face) in direction of cell centre
// so it is clearly inside.
Foam::point Foam::sampledSet::pushIn
(
    const point& facePt,
    const label faceI
) const
{
    label cellI = mesh().faceOwner()[faceI];
    const point& cC = mesh().cellCentres()[cellI];

    point newPosition = facePt;

    // Taken from particle::initCellFacePt()
    label tetFaceI;
    label tetPtI;
    mesh().findTetFacePt(cellI, facePt, tetFaceI, tetPtI);

    if (tetFaceI == -1 || tetPtI == -1)
    {
        newPosition = facePt;

        label trap(1.0/particle::trackingCorrectionTol + 1);

        label iterNo = 0;

        do
        {
            newPosition += particle::trackingCorrectionTol*(cC - facePt);

            mesh().findTetFacePt
            (
                cellI,
                newPosition,
                tetFaceI,
                tetPtI
            );

            iterNo++;

        } while (tetFaceI < 0  && iterNo <= trap);
    }

    if (tetFaceI == -1)
    {
        FatalErrorIn
        (
            "sampledSet::pushIn(const point&, const label)"
        )   << "After pushing " << facePt << " to " << newPosition
            << " it is still outside face " << faceI
            << " at " << mesh().faceCentres()[faceI]
            << " of cell " << cellI
            << " at " << cC << endl
            << "Please change your starting point"
            << abort(FatalError);
    }

    //Info<< "pushIn : moved " << facePt << " to " << newPosition
    //    << endl;

    return newPosition;
}


// Calculates start of tracking given samplePt and first boundary intersection
// (bPoint, bFaceI). bFaceI == -1 if no boundary intersection.
// Returns true if trackPt is sampling point
bool Foam::sampledSet::getTrackingPoint
(
    const vector& offset,
    const point& samplePt,
    const point& bPoint,
    const label bFaceI,

    point& trackPt,
    label& trackCellI,
    label& trackFaceI
) const
{
    const scalar smallDist = mag(tol*offset);

    bool isGoodSample = false;

    if (bFaceI == -1)
    {
        // No boundary intersection. Try and find cell samplePt is in
        trackCellI = mesh().findCell(samplePt, searchEngine_.decompMode());

        if
        (
            (trackCellI == -1)
        || !mesh().pointInCell
            (
                samplePt,
                trackCellI,
                searchEngine_.decompMode()
            )
        )
        {
            // Line samplePt - end_ does not intersect domain at all.
            // (or is along edge)
            //Info<< "getTrackingPoint : samplePt outside domain : "
            //    << "  samplePt:" << samplePt
            //    << endl;

            trackCellI = -1;
            trackFaceI = -1;

            isGoodSample = false;
        }
        else
        {
            // start is inside. Use it as tracking point
            //Info<< "getTrackingPoint : samplePt inside :"
            //    << "  samplePt:" << samplePt
            //    << "  trackCellI:" << trackCellI
            //    << endl;

            trackPt = samplePt;
            trackFaceI = -1;

            isGoodSample = true;
        }
    }
    else if (mag(samplePt - bPoint) < smallDist)
    {
        //Info<< "getTrackingPoint : samplePt:" << samplePt
        //    << " close to bPoint:"
        //    << bPoint << endl;

        // samplePt close to bPoint. Snap to it
        trackPt = pushIn(bPoint, bFaceI);
        trackFaceI = bFaceI;
        trackCellI = getBoundaryCell(trackFaceI);

        isGoodSample = true;
    }
    else
    {
        scalar sign = calcSign(bFaceI, samplePt);

        if (sign < 0)
        {
            // samplePt inside or marginally outside.
            trackPt = samplePt;
            trackFaceI = -1;
            trackCellI = mesh().findCell(trackPt, searchEngine_.decompMode());

            isGoodSample = true;
        }
        else
        {
            // samplePt outside. use bPoint
            trackPt = pushIn(bPoint, bFaceI);
            trackFaceI = bFaceI;
            trackCellI = getBoundaryCell(trackFaceI);

            isGoodSample = false;
        }
    }

    if (debug)
    {
        Info<< "sampledSet::getTrackingPoint :"
            << " offset:" << offset
            << " samplePt:" << samplePt
            << " bPoint:" << bPoint
            << " bFaceI:" << bFaceI
            << endl << "   Calculated first tracking point :"
            << " trackPt:" << trackPt
            << " trackCellI:" << trackCellI
            << " trackFaceI:" << trackFaceI
            << " isGoodSample:" << isGoodSample
            << endl;
    }

    return isGoodSample;
}


void Foam::sampledSet::setSamples
(
    const List<point>& samplingPts,
    const labelList& samplingCells,
    const labelList& samplingFaces,
    const labelList& samplingSegments,
    const scalarList& samplingCurveDist
)
{
    setSize(samplingPts.size());
    cells_.setSize(samplingCells.size());
    faces_.setSize(samplingFaces.size());
    segments_.setSize(samplingSegments.size());
    curveDist_.setSize(samplingCurveDist.size());

    if
    (
        (cells_.size() != size())
     || (faces_.size() != size())
     || (segments_.size() != size())
     || (curveDist_.size() != size())
    )
    {
        FatalErrorIn("sampledSet::setSamples()")
            << "sizes not equal : "
            << "  points:" << size()
            << "  cells:" << cells_.size()
            << "  faces:" << faces_.size()
            << "  segments:" << segments_.size()
            << "  curveDist:" << curveDist_.size()
            << abort(FatalError);
    }

    forAll(samplingPts, sampleI)
    {
        operator[](sampleI) = samplingPts[sampleI];
    }
    curveDist_ = samplingCurveDist;

    cells_ = samplingCells;
    faces_ = samplingFaces;
    segments_ = samplingSegments;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSet::sampledSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis
)
:
    coordSet(name, axis),
    mesh_(mesh),
    searchEngine_(searchEngine),
    segments_(0),
    cells_(0),
    faces_(0)
{}


Foam::sampledSet::sampledSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    coordSet(name, dict.lookup("axis")),
    mesh_(mesh),
    searchEngine_(searchEngine),
    segments_(0),
    cells_(0),
    faces_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSet::~sampledSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::sampledSet> Foam::sampledSet::New
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
{
    const word sampleType(dict.lookup("type"));

    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(sampleType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "sampledSet::New"
            "(const word&, const polyMesh&, const meshSearch&"
            ", const dictionary&)"
        )   << "Unknown sample type "
            << sampleType << nl << nl
            << "Valid sample types : " << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<sampledSet>
    (
        cstrIter()
        (
            name,
            mesh,
            searchEngine,
            dict
        )
    );
}


Foam::Ostream& Foam::sampledSet::write(Ostream& os) const
{
    coordSet::write(os);

    os  << endl << "\t(cellI)\t(faceI)" << endl;

    forAll(*this, sampleI)
    {
        os  << '\t' << cells_[sampleI]
            << '\t' << faces_[sampleI]
            << endl;
    }

    return os;
}


// ************************************************************************* //
