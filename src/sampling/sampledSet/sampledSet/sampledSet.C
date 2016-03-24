/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    defineTypeNameAndDebug(sampledSet, 0);
    defineRunTimeSelectionTable(sampledSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::sampledSet::getBoundaryCell(const label faceI) const
{
    return mesh().faceOwner()[faceI];
}


Foam::label Foam::sampledSet::getNeighbourCell(const label faceI) const
{
    if (faceI >= mesh().nInternalFaces())
    {
        return mesh().faceOwner()[faceI];
    }
    else
    {
        return mesh().faceNeighbour()[faceI];
    }
}


Foam::label Foam::sampledSet::pointInCell
(
    const point& p,
    const label samplei
) const
{
    // Collect the face owner and neighbour cells of the sample into an array
    // for convenience
    label cells[4] =
    {
        mesh().faceOwner()[faces_[samplei]],
        getNeighbourCell(faces_[samplei]),
        mesh().faceOwner()[faces_[samplei+1]],
        getNeighbourCell(faces_[samplei+1])
    };

    // Find the sampled cell by checking the owners and neighbours of the
    // sampled faces
    label cellm =
        (cells[0] == cells[2] || cells[0] == cells[3]) ? cells[0]
      : (cells[1] == cells[2] || cells[1] == cells[3]) ? cells[1]
      : -1;

    if (cellm != -1)
    {
        // If found the sampled cell check the point is in the cell
        // otherwise ignore
        if (!mesh().pointInCell(p, cellm, searchEngine_.decompMode()))
        {
           cellm = -1;

            if (debug)
            {
                WarningInFunction
                    << "Could not find mid-point " << p
                    << " cell " << cellm << endl;
            }
        }
    }
    else
    {
        // If the sample does not pass through a single cell check if the point
        // is in any of the owners or neighbours otherwise ignore
        for (label i=0; i<4; i++)
        {
            if (mesh().pointInCell(p, cells[i], searchEngine_.decompMode()))
            {
                return cells[i];
            }
        }

        if (debug)
        {
            WarningInFunction
                << "Could not find cell for mid-point" << nl
                << "  samplei: " << samplei
                << "  pts[samplei]: " << operator[](samplei)
                << "  face[samplei]: " << faces_[samplei]
                << "  pts[samplei+1]: " << operator[](samplei+1)
                << "  face[samplei+1]: " << faces_[samplei+1]
                << "  cellio: " << cells[0]
                << "  cellin: " << cells[1]
                << "  celljo: " << cells[2]
                << "  celljn: " << cells[3]
                << endl;
        }
    }

    return cellm;
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
        // Sample on face centre. Regard as inside
        return -1;
    }

    vec /= magVec;

    vector n = mesh().faceAreas()[faceI];

    n /= mag(n) + VSMALL;

    return n & vec;
}


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
        FatalErrorInFunction
            << "After pushing " << facePt << " to " << newPosition
            << " it is still outside face " << faceI
            << " at " << mesh().faceCentres()[faceI]
            << " of cell " << cellI
            << " at " << cC << endl
            << "Please change your starting point"
            << abort(FatalError);
    }

    return newPosition;
}


bool Foam::sampledSet::getTrackingPoint
(
    const point& samplePt,
    const point& bPoint,
    const label bFaceI,
    const scalar smallDist,

    point& trackPt,
    label& trackCellI,
    label& trackFaceI
) const
{
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

            trackCellI = -1;
            trackFaceI = -1;

            isGoodSample = false;
        }
        else
        {
            // Start is inside. Use it as tracking point

            trackPt = samplePt;
            trackFaceI = -1;

            isGoodSample = true;
        }
    }
    else if (mag(samplePt - bPoint) < smallDist)
    {
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
        InfoInFunction
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
        FatalErrorInFunction
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
        FatalErrorInFunction
            << "Unknown sample type "
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
