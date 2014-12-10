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

#include "searchablePlane.H"
#include "addToRunTimeSelectionTable.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(searchablePlane, 0);
addToRunTimeSelectionTable(searchableSurface, searchablePlane, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::pointIndexHit Foam::searchablePlane::findLine
(
    const point& start,
    const point& end
) const
{
    pointIndexHit info(true, vector::zero, 0);

    linePointRef l(start, end);

    scalar t = lineIntersect(l);

    if (t < 0 || t > 1)
    {
        info.setMiss();
        info.setIndex(-1);
    }
    else
    {
        info.setPoint(start+t*l.vec());
    }

    return info;
}


Foam::boundBox Foam::searchablePlane::calcBounds() const
{
    point max(VGREAT, VGREAT, VGREAT);

    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (mag(normal()[dir]) - 1 < SMALL)
        {
            max[dir] = 0;

            break;
        }
    }

    point min = -max;

    return boundBox(min, max);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchablePlane::searchablePlane
(
    const IOobject& io,
    const point& basePoint,
    const vector& normal
)
:
    searchableSurface(io),
    plane(basePoint, normal)
{
    bounds() = calcBounds();
}


Foam::searchablePlane::searchablePlane
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    plane(dict)
{
    bounds() = calcBounds();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchablePlane::~searchablePlane()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchablePlane::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(1);
        regions_[0] = "region0";
    }
    return regions_;
}


void Foam::searchablePlane::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres.setSize(1);
    centres[0] = refPoint();

    radiusSqr.setSize(1);
    radiusSqr[0] = Foam::sqr(GREAT);
}


void Foam::searchablePlane::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    info.setSize(samples.size());

    forAll(samples, i)
    {
        info[i].setPoint(nearestPoint(samples[i]));

        if (magSqr(samples[i]-info[i].rawPoint()) > nearestDistSqr[i])
        {
            info[i].setIndex(-1);
            info[i].setMiss();
        }
        else
        {
            info[i].setIndex(0);
            info[i].setHit();
        }
    }
}


void Foam::searchablePlane::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        info[i] = findLine(start[i], end[i]);
    }
}


void Foam::searchablePlane::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    findLine(start, end, info);
}


void Foam::searchablePlane::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit> >& info
) const
{
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


void Foam::searchablePlane::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.setSize(info.size());
    region = 0;
}


void Foam::searchablePlane::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& n
) const
{
    n.setSize(info.size());
    n = normal();
}


void Foam::searchablePlane::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    FatalErrorIn
    (
        "searchableCollection::getVolumeType(const pointField&"
        ", List<volumeType>&) const"
    )   << "Volume type not supported for plane."
        << exit(FatalError);
}


// ************************************************************************* //
