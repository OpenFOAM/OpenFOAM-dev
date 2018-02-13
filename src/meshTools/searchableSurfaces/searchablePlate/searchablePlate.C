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

#include "searchablePlate.H"
#include "addToRunTimeSelectionTable.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchablePlate, 0);
    addToRunTimeSelectionTable(searchableSurface, searchablePlate, dict);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::direction Foam::searchablePlate::calcNormal(const point& span)
{
    direction normalDir = 3;

    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (span[dir] < 0)
        {
            FatalErrorInFunction
                << "Span should have two positive and one zero entry. Now:"
                << span << exit(FatalError);
        }
        else if (span[dir] < vSmall)
        {
            if (normalDir == 3)
            {
                normalDir = dir;
            }
            else
            {
                // Multiple zero entries. Flag and exit.
                normalDir = 3;
                break;
            }
        }
    }

    if (normalDir == 3)
    {
        FatalErrorInFunction
            << "Span should have two positive and one zero entry. Now:"
            << span << exit(FatalError);
    }

    return normalDir;
}


// Returns miss or hit with face (always 0)
Foam::pointIndexHit Foam::searchablePlate::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    // For every component direction can be
    // left of min, right of max or in between.
    // - outside points: project first one x plane (either min().x()
    // or max().x()), then onto y plane and finally z. You should be left
    // with intersection point
    // - inside point: find nearest side (compare to mid point). Project onto
    //   that.

    // Project point on plane.
    pointIndexHit info(true, sample, 0);
    info.rawPoint()[normalDir_] = origin_[normalDir_];

    // Clip to edges if outside
    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (dir != normalDir_)
        {
            if (info.rawPoint()[dir] < origin_[dir])
            {
                info.rawPoint()[dir] = origin_[dir];
            }
            else if (info.rawPoint()[dir] > origin_[dir]+span_[dir])
            {
                info.rawPoint()[dir] = origin_[dir]+span_[dir];
            }
        }
    }

    // Check if outside. Optimisation: could do some checks on distance already
    // on components above
    if (magSqr(info.rawPoint() - sample) > nearestDistSqr)
    {
        info.setMiss();
        info.setIndex(-1);
    }

    return info;
}


Foam::pointIndexHit Foam::searchablePlate::findLine
(
    const point& start,
    const point& end
) const
{
    pointIndexHit info
    (
        true,
        Zero,
        0
    );

    const vector dir(end-start);

    if (mag(dir[normalDir_]) < vSmall)
    {
        info.setMiss();
        info.setIndex(-1);
    }
    else
    {
        scalar t = (origin_[normalDir_]-start[normalDir_]) / dir[normalDir_];

        if (t < 0 || t > 1)
        {
            info.setMiss();
            info.setIndex(-1);
        }
        else
        {
            info.rawPoint() = start+t*dir;
            info.rawPoint()[normalDir_] = origin_[normalDir_];

            // Clip to edges
            for (direction dir = 0; dir < vector::nComponents; dir++)
            {
                if (dir != normalDir_)
                {
                    if (info.rawPoint()[dir] < origin_[dir])
                    {
                        info.setMiss();
                        info.setIndex(-1);
                        break;
                    }
                    else if (info.rawPoint()[dir] > origin_[dir]+span_[dir])
                    {
                        info.setMiss();
                        info.setIndex(-1);
                        break;
                    }
                }
            }
        }
    }

    // Debug
    if (info.hit())
    {
        treeBoundBox bb(origin_, origin_+span_);
        bb.min()[normalDir_] -= 1e-6;
        bb.max()[normalDir_] += 1e-6;

        if (!bb.contains(info.hitPoint()))
        {
            FatalErrorInFunction
                << "bb:" << bb << endl
                << "origin_:" << origin_ << endl
                << "span_:" << span_ << endl
                << "normalDir_:" << normalDir_ << endl
                << "hitPoint:" << info.hitPoint()
                << abort(FatalError);
        }
    }

    return info;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchablePlate::searchablePlate
(
    const IOobject& io,
    const point& origin,
    const vector& span
)
:
    searchableSurface(io),
    origin_(origin),
    span_(span),
    normalDir_(calcNormal(span_))
{
    if (debug)
    {
        InfoInFunction
            << " origin:" << origin_
            << " origin+span:" << origin_+span_
            << " normal:" << vector::componentNames[normalDir_]
            << endl;
    }

    bounds() = boundBox(origin_, origin_ + span_);
}


Foam::searchablePlate::searchablePlate
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    origin_(dict.lookup("origin")),
    span_(dict.lookup("span")),
    normalDir_(calcNormal(span_))
{
    if (debug)
    {
        InfoInFunction
            << " origin:" << origin_
            << " origin+span:" << origin_+span_
            << " normal:" << vector::componentNames[normalDir_]
            << endl;
    }

    bounds() = boundBox(origin_, origin_ + span_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchablePlate::~searchablePlate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchablePlate::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(1);
        regions_[0] = "region0";
    }
    return regions_;
}


Foam::tmp<Foam::pointField> Foam::searchablePlate::coordinates() const
{
    return tmp<pointField>(new pointField(1, origin_ + 0.5*span_));
}


void Foam::searchablePlate::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres.setSize(1);
    centres[0] = origin_ + 0.5*span_;

    radiusSqr.setSize(1);
    radiusSqr[0] = Foam::magSqr(0.5*span_);

    // Add a bit to make sure all points are tested inside
    radiusSqr += Foam::sqr(small);
}


Foam::tmp<Foam::pointField> Foam::searchablePlate::points() const
{
    tmp<pointField> tPts(new pointField(4));
    pointField& pts = tPts.ref();

    pts[0] = origin_;
    pts[2] = origin_ + span_;

    if (span_.x() < span_.y() && span_.x() < span_.z())
    {
        pts[1] = origin_ + point(0, span_.y(), 0);
        pts[3] = origin_ + point(0, 0, span_.z());
    }
    else if (span_.y() < span_.z())
    {
        pts[1] = origin_ + point(span_.x(), 0, 0);
        pts[3] = origin_ + point(0, 0, span_.z());
    }
    else
    {
        pts[1] = origin_ + point(span_.x(), 0, 0);
        pts[3] = origin_ + point(0, span_.y(), 0);
    }

    return tPts;
}


bool Foam::searchablePlate::overlaps(const boundBox& bb) const
{
    return
    (
        (origin_.x() + span_.x()) >= bb.min().x()
     && origin_.x() <= bb.max().x()
     && (origin_.y() + span_.y()) >= bb.min().y()
     && origin_.y() <= bb.max().y()
     && (origin_.z() + span_.z()) >= bb.min().z()
     && origin_.z() <= bb.max().z()
    );
}


void Foam::searchablePlate::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    info.setSize(samples.size());

    forAll(samples, i)
    {
        info[i] = findNearest(samples[i], nearestDistSqr[i]);
    }
}


void Foam::searchablePlate::findLine
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


void Foam::searchablePlate::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    findLine(start, end, info);
}


void Foam::searchablePlate::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    List<pointIndexHit> nearestInfo;
    findLine(start, end, nearestInfo);

    info.setSize(start.size());
    forAll(info, pointi)
    {
        if (nearestInfo[pointi].hit())
        {
            info[pointi].setSize(1);
            info[pointi][0] = nearestInfo[pointi];
        }
        else
        {
            info[pointi].clear();
        }
    }
}


void Foam::searchablePlate::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.setSize(info.size());
    region = 0;
}


void Foam::searchablePlate::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    normal.setSize(info.size());
    normal = Zero;
    forAll(normal, i)
    {
        normal[i][normalDir_] = 1.0;
    }
}


void Foam::searchablePlate::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    FatalErrorInFunction
        << "Volume type not supported for plate."
        << exit(FatalError);
}


// ************************************************************************* //
