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

#include "boundBox.H"
#include "PstreamReduceOps.H"
#include "tmp.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::boundBox::great(vGreat);

const Foam::boundBox Foam::boundBox::greatBox
(
    point(-vGreat, -vGreat, -vGreat),
    point(vGreat, vGreat, vGreat)
);


const Foam::boundBox Foam::boundBox::invertedBox
(
    point(vGreat, vGreat, vGreat),
    point(-vGreat, -vGreat, -vGreat)
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::boundBox::calculate(const UList<point>& points, const bool doReduce)
{
    if (points.empty())
    {
        min_ = Zero;
        max_ = Zero;

        if (doReduce && Pstream::parRun())
        {
            // Use values that get overwritten by reduce minOp, maxOp below
            min_ = point(vGreat, vGreat, vGreat);
            max_ = point(-vGreat, -vGreat, -vGreat);
        }
    }
    else
    {
        min_ = points[0];
        max_ = points[0];


        for (label i = 1; i < points.size(); i++)
        {
            min_ = ::Foam::min(min_, points[i]);
            max_ = ::Foam::max(max_, points[i]);
        }
    }

    // Reduce parallel information
    if (doReduce)
    {
        reduce(min_, minOp<point>());
        reduce(max_, maxOp<point>());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundBox::boundBox(const UList<point>& points, const bool doReduce)
:
    min_(Zero),
    max_(Zero)
{
    calculate(points, doReduce);
}


Foam::boundBox::boundBox(const tmp<pointField>& points, const bool doReduce)
:
    min_(Zero),
    max_(Zero)
{
    calculate(points(), doReduce);
    points.clear();
}


Foam::boundBox::boundBox
(
    const UList<point>& points,
    const labelUList& indices,
    const bool doReduce
)
:
    min_(Zero),
    max_(Zero)
{
    if (points.empty() || indices.empty())
    {
        if (doReduce && Pstream::parRun())
        {
            // Use values that get overwritten by reduce minOp, maxOp below
            min_ = point(vGreat, vGreat, vGreat);
            max_ = point(-vGreat, -vGreat, -vGreat);
        }
    }
    else
    {
        min_ = points[indices[0]];
        max_ = points[indices[0]];

        for (label i=1; i < indices.size(); ++i)
        {
            min_ = ::Foam::min(min_, points[indices[i]]);
            max_ = ::Foam::max(max_, points[indices[i]]);
        }
    }

    // Reduce parallel information
    if (doReduce)
    {
        reduce(min_, minOp<point>());
        reduce(max_, maxOp<point>());
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::boundBox::points() const
{
    tmp<pointField> tPts = tmp<pointField>(new pointField(8));
    pointField& pt = tPts.ref();

    pt[0] = min_;                                   // min-x, min-y, min-z
    pt[1] = point(max_.x(), min_.y(), min_.z());    // max-x, min-y, min-z
    pt[2] = point(max_.x(), max_.y(), min_.z());    // max-x, max-y, min-z
    pt[3] = point(min_.x(), max_.y(), min_.z());    // min-x, max-y, min-z
    pt[4] = point(min_.x(), min_.y(), max_.z());    // min-x, min-y, max-z
    pt[5] = point(max_.x(), min_.y(), max_.z());    // max-x, min-y, max-z
    pt[6] = max_;                                   // max-x, max-y, max-z
    pt[7] = point(min_.x(), max_.y(), max_.z());    // min-x, max-y, max-z

    return tPts;
}


Foam::faceList Foam::boundBox::faces()
{
    faceList faces(6);

    forAll(faces, fI)
    {
        faces[fI].setSize(4);
    }

    faces[0][0] = 0;
    faces[0][1] = 1;
    faces[0][2] = 2;
    faces[0][3] = 3;

    faces[1][0] = 2;
    faces[1][1] = 6;
    faces[1][2] = 7;
    faces[1][3] = 3;

    faces[2][0] = 0;
    faces[2][1] = 4;
    faces[2][2] = 5;
    faces[2][3] = 1;

    faces[3][0] = 4;
    faces[3][1] = 7;
    faces[3][2] = 6;
    faces[3][3] = 5;

    faces[4][0] = 3;
    faces[4][1] = 7;
    faces[4][2] = 4;
    faces[4][3] = 0;

    faces[5][0] = 1;
    faces[5][1] = 5;
    faces[5][2] = 6;
    faces[5][3] = 2;

    return faces;
}


void Foam::boundBox::inflate(const scalar s)
{
    vector ext = vector::one*s*mag();

    min_ -= ext;
    max_ += ext;
}


bool Foam::boundBox::contains(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    forAll(points, i)
    {
        if (!contains(points[i]))
        {
            return false;
        }
    }

    return true;
}


bool Foam::boundBox::contains
(
    const UList<point>& points,
    const labelUList& indices
) const
{
    if (points.empty() || indices.empty())
    {
        return true;
    }

    forAll(indices, i)
    {
        if (!contains(points[indices[i]]))
        {
            return false;
        }
    }

    return true;
}


bool Foam::boundBox::containsAny(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    forAll(points, i)
    {
        if (contains(points[i]))
        {
            return true;
        }
    }

    return false;
}


bool Foam::boundBox::containsAny
(
    const UList<point>& points,
    const labelUList& indices
) const
{
    if (points.empty() || indices.empty())
    {
        return true;
    }

    forAll(indices, i)
    {
        if (contains(points[indices[i]]))
        {
            return true;
        }
    }

    return false;
}


Foam::point Foam::boundBox::nearest(const point& pt) const
{
    // Clip the point to the range of the bounding box
    const scalar surfPtx = Foam::max(Foam::min(pt.x(), max_.x()), min_.x());
    const scalar surfPty = Foam::max(Foam::min(pt.y(), max_.y()), min_.y());
    const scalar surfPtz = Foam::max(Foam::min(pt.z(), max_.z()), min_.z());

    return point(surfPtx, surfPty, surfPtz);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const boundBox& bb)
{
    if (os.format() == IOstream::ASCII)
    {
        os << bb.min_ << token::SPACE << bb.max_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&bb.min_),
            sizeof(boundBox)
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const boundBox&)");
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, boundBox& bb)
{
    if (is.format() == IOstream::ASCII)
    {
        is >> bb.min_ >> bb.max_;
    }
    else
    {
        is.read
        (
            reinterpret_cast<char*>(&bb.min_),
            sizeof(boundBox)
        );
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, boundBox&)");
    return is;
}


// ************************************************************************* //
