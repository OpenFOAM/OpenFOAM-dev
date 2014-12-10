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

#include "treeBoundBox.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::treeBoundBox::great(GREAT);

const Foam::treeBoundBox Foam::treeBoundBox::greatBox
(
    vector(-GREAT, -GREAT, -GREAT),
    vector(GREAT, GREAT, GREAT)
);


const Foam::treeBoundBox Foam::treeBoundBox::invertedBox
(
    vector(GREAT, GREAT, GREAT),
    vector(-GREAT, -GREAT, -GREAT)
);


//! \cond - skip documentation : local scope only
const Foam::label facesArray[6][4] =
{
    {0, 4, 6, 2}, // left
    {1, 3, 7, 5}, // right
    {0, 1, 5, 4}, // bottom
    {2, 6, 7, 3}, // top
    {0, 2, 3, 1}, // back
    {4, 5, 7, 6}  // front
};
//! \endcond


const Foam::faceList Foam::treeBoundBox::faces
(
    initListList<face, label, 6, 4>(facesArray)
);


//! \cond - skip documentation : local scope only
const Foam::label edgesArray[12][2] =
{
    {0, 1}, // 0
    {1, 3},
    {2, 3}, // 2
    {0, 2},
    {4, 5}, // 4
    {5, 7},
    {6, 7}, // 6
    {4, 6},
    {0, 4}, // 8
    {1, 5},
    {3, 7}, // 10
    {2, 6}
};
//! \endcond


const Foam::edgeList Foam::treeBoundBox::edges
(
    //initListList<edge, label, 12, 2>(edgesArray)
    calcEdges(edgesArray)
);


const Foam::FixedList<Foam::vector, 6> Foam::treeBoundBox::faceNormals
(
    calcFaceNormals()
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::edgeList Foam::treeBoundBox::calcEdges(const label edgesArray[12][2])
{
    edgeList edges(12);
    forAll(edges, edgeI)
    {
        edges[edgeI][0] = edgesArray[edgeI][0];
        edges[edgeI][1] = edgesArray[edgeI][1];
    }
    return edges;
}


Foam::FixedList<Foam::vector, 6> Foam::treeBoundBox::calcFaceNormals()
{
    FixedList<vector, 6> normals;
    normals[LEFT]   = vector(-1,  0,  0);
    normals[RIGHT]  = vector( 1,  0,  0);
    normals[BOTTOM] = vector( 0, -1,  0);
    normals[TOP]    = vector( 0,  1,  0);
    normals[BACK]   = vector( 0,  0, -1);
    normals[FRONT]  = vector( 0,  0,  1);
    return normals;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::treeBoundBox::treeBoundBox(const UList<point>& points)
:
    boundBox(points, false)
{
    if (points.empty())
    {
        WarningIn
        (
            "treeBoundBox::treeBoundBox(const UList<point>&)"
        )   << "cannot find bounding box for zero-sized pointField, "
            << "returning zero" << endl;

        return;
    }
}


Foam::treeBoundBox::treeBoundBox
(
    const UList<point>& points,
    const labelUList& indices
)
:
    boundBox(points, indices, false)
{
    if (points.empty() || indices.empty())
    {
        WarningIn
        (
            "treeBoundBox::treeBoundBox"
            "(const UList<point>&, const labelUList&)"
        )   << "cannot find bounding box for zero-sized pointField, "
            << "returning zero" << endl;

        return;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::treeBoundBox::points() const
{
    tmp<pointField> tPts = tmp<pointField>(new pointField(8));

    pointField& points = tPts();

    forAll(points, octant)
    {
        points[octant] = corner(octant);
    }

    return tPts;
}


Foam::treeBoundBox Foam::treeBoundBox::subBbox(const direction octant) const
{
    return subBbox(midpoint(), octant);
}


// Octant to bounding box using permutation only.
Foam::treeBoundBox Foam::treeBoundBox::subBbox
(
    const point& mid,
    const direction octant
) const
{
    if (octant > 7)
    {
        FatalErrorIn
        (
            "treeBoundBox::subBbox(const point&, const direction)"
        )   << "octant should be [0..7]"
            << abort(FatalError);
    }

    // start with a copy of this bounding box and adjust limits accordingly
    treeBoundBox subBb(*this);
    point& bbMin = subBb.min();
    point& bbMax = subBb.max();

    if (octant & treeBoundBox::RIGHTHALF)
    {
        bbMin.x() = mid.x();    // mid -> max
    }
    else
    {
        bbMax.x() = mid.x();    // min -> mid
    }

    if (octant & treeBoundBox::TOPHALF)
    {
        bbMin.y() = mid.y();    // mid -> max
    }
    else
    {
        bbMax.y() = mid.y();    // min -> mid
    }

    if (octant & treeBoundBox::FRONTHALF)
    {
        bbMin.z() = mid.z();    // mid -> max
    }
    else
    {
        bbMax.z() = mid.z();    // min -> mid
    }

    return subBb;
}


// line intersection. Returns true if line (start to end) inside
// bb or intersects bb. Sets pt to intersection.
//
// Sutherlands algorithm:
//   loop
//     - start = intersection of line with one of the planes bounding
//       the bounding box
//     - stop if start inside bb (return true)
//     - stop if start and end in same 'half' (e.g. both above bb)
//       (return false)
//
// Uses posBits to efficiently determine 'half' in which start and end
// point are.
//
// Note:
//   - sets coordinate to exact position: e.g. pt.x() = min().x()
//     since plane intersect routine might have truncation error.
//     This makes sure that posBits tests 'inside'
bool Foam::treeBoundBox::intersects
(
    const point& overallStart,
    const vector& overallVec,
    const point& start,
    const point& end,
    point& pt,
    direction& ptOnFaces
) const
{
    const direction endBits = posBits(end);
    pt = start;

    // Allow maximum of 3 clips.
    for (label i = 0; i < 4; ++i)
    {
        direction ptBits = posBits(pt);

        if (ptBits == 0)
        {
            // pt inside bb
            ptOnFaces = faceBits(pt);
            return true;
        }

        if ((ptBits & endBits) != 0)
        {
            // pt and end in same block outside of bb
            ptOnFaces = faceBits(pt);
            return false;
        }

        if (ptBits & LEFTBIT)
        {
            // Intersect with plane V=min, n=-1,0,0
            if (Foam::mag(overallVec.x()) > VSMALL)
            {
                scalar s = (min().x() - overallStart.x())/overallVec.x();
                pt.x() = min().x();
                pt.y() = overallStart.y() + overallVec.y()*s;
                pt.z() = overallStart.z() + overallVec.z()*s;
            }
            else
            {
                // Vector not in x-direction. But still intersecting bb planes.
                // So must be close - just snap to bb.
                pt.x() = min().x();
            }
        }
        else if (ptBits & RIGHTBIT)
        {
            // Intersect with plane V=max, n=1,0,0
            if (Foam::mag(overallVec.x()) > VSMALL)
            {
                scalar s = (max().x() - overallStart.x())/overallVec.x();
                pt.x() = max().x();
                pt.y() = overallStart.y() + overallVec.y()*s;
                pt.z() = overallStart.z() + overallVec.z()*s;
            }
            else
            {
                pt.x() = max().x();
            }
        }
        else if (ptBits & BOTTOMBIT)
        {
            // Intersect with plane V=min, n=0,-1,0
            if (Foam::mag(overallVec.y()) > VSMALL)
            {
                scalar s = (min().y() - overallStart.y())/overallVec.y();
                pt.x() = overallStart.x() + overallVec.x()*s;
                pt.y() = min().y();
                pt.z() = overallStart.z() + overallVec.z()*s;
            }
            else
            {
                pt.x() = min().y();
            }
        }
        else if (ptBits & TOPBIT)
        {
            // Intersect with plane V=max, n=0,1,0
            if (Foam::mag(overallVec.y()) > VSMALL)
            {
                scalar s = (max().y() - overallStart.y())/overallVec.y();
                pt.x() = overallStart.x() + overallVec.x()*s;
                pt.y() = max().y();
                pt.z() = overallStart.z() + overallVec.z()*s;
            }
            else
            {
                pt.y() = max().y();
            }
        }
        else if (ptBits & BACKBIT)
        {
            // Intersect with plane V=min, n=0,0,-1
            if (Foam::mag(overallVec.z()) > VSMALL)
            {
                scalar s = (min().z() - overallStart.z())/overallVec.z();
                pt.x() = overallStart.x() + overallVec.x()*s;
                pt.y() = overallStart.y() + overallVec.y()*s;
                pt.z() = min().z();
            }
            else
            {
                pt.z() = min().z();
            }
        }
        else if (ptBits & FRONTBIT)
        {
            // Intersect with plane V=max, n=0,0,1
            if (Foam::mag(overallVec.z()) > VSMALL)
            {
                scalar s = (max().z() - overallStart.z())/overallVec.z();
                pt.x() = overallStart.x() + overallVec.x()*s;
                pt.y() = overallStart.y() + overallVec.y()*s;
                pt.z() = max().z();
            }
            else
            {
                pt.z() = max().z();
            }
        }
    }

    // Can end up here if the end point is on the edge of the boundBox
    return true;
}


bool Foam::treeBoundBox::intersects
(
    const point& start,
    const point& end,
    point& pt
) const
{
    direction ptBits;
    return intersects(start, end-start, start, end, pt, ptBits);
}


bool Foam::treeBoundBox::contains(const vector& dir, const point& pt) const
{
    //
    // Compare all components against min and max of bb
    //
    for (direction cmpt=0; cmpt<3; cmpt++)
    {
        if (pt[cmpt] < min()[cmpt])
        {
            return false;
        }
        else if (pt[cmpt] == min()[cmpt])
        {
            // On edge. Outside if direction points outwards.
            if (dir[cmpt] < 0)
            {
                return false;
            }
        }

        if (pt[cmpt] > max()[cmpt])
        {
            return false;
        }
        else if (pt[cmpt] == max()[cmpt])
        {
            // On edge. Outside if direction points outwards.
            if (dir[cmpt] > 0)
            {
                return false;
            }
        }
    }

    // All components inside bb
    return true;
}


// Code position of pt on bounding box faces
Foam::direction Foam::treeBoundBox::faceBits(const point& pt) const
{
    direction faceBits = 0;
    if (pt.x() == min().x())
    {
        faceBits |= LEFTBIT;
    }
    else if (pt.x() == max().x())
    {
        faceBits |= RIGHTBIT;
    }

    if (pt.y() == min().y())
    {
        faceBits |= BOTTOMBIT;
    }
    else if (pt.y() == max().y())
    {
        faceBits |= TOPBIT;
    }

    if (pt.z() == min().z())
    {
        faceBits |= BACKBIT;
    }
    else if (pt.z() == max().z())
    {
        faceBits |= FRONTBIT;
    }
    return faceBits;
}


// Code position of point relative to box
Foam::direction Foam::treeBoundBox::posBits(const point& pt) const
{
    direction posBits = 0;

    if (pt.x() < min().x())
    {
        posBits |= LEFTBIT;
    }
    else if (pt.x() > max().x())
    {
        posBits |= RIGHTBIT;
    }

    if (pt.y() < min().y())
    {
        posBits |= BOTTOMBIT;
    }
    else if (pt.y() > max().y())
    {
        posBits |= TOPBIT;
    }

    if (pt.z() < min().z())
    {
        posBits |= BACKBIT;
    }
    else if (pt.z() > max().z())
    {
        posBits |= FRONTBIT;
    }
    return posBits;
}


// nearest and furthest corner coordinate.
// !names of treeBoundBox::min() and treeBoundBox::max() are confusing!
void Foam::treeBoundBox::calcExtremities
(
    const point& pt,
    point& nearest,
    point& furthest
) const
{
    scalar nearX, nearY, nearZ;
    scalar farX, farY, farZ;

    if (Foam::mag(min().x() - pt.x()) < Foam::mag(max().x() - pt.x()))
    {
        nearX = min().x();
        farX = max().x();
    }
    else
    {
        nearX = max().x();
        farX = min().x();
    }

    if (Foam::mag(min().y() - pt.y()) < Foam::mag(max().y() - pt.y()))
    {
        nearY = min().y();
        farY = max().y();
    }
    else
    {
        nearY = max().y();
        farY = min().y();
    }

    if (Foam::mag(min().z() - pt.z()) < Foam::mag(max().z() - pt.z()))
    {
        nearZ = min().z();
        farZ = max().z();
    }
    else
    {
        nearZ = max().z();
        farZ = min().z();
    }

    nearest = point(nearX, nearY, nearZ);
    furthest = point(farX, farY, farZ);
}


Foam::scalar Foam::treeBoundBox::maxDist(const point& pt) const
{
    point near, far;
    calcExtremities(pt, near, far);

    return Foam::mag(far - pt);
}


// Distance comparator
// Compare all vertices of bounding box against all of other bounding
// box to see if all vertices of one are nearer
Foam::label Foam::treeBoundBox::distanceCmp
(
    const point& pt,
    const treeBoundBox& other
) const
{
    //
    // Distance point <-> nearest and furthest away vertex of this
    //

    point nearThis, farThis;

    // get nearest and furthest away vertex
    calcExtremities(pt, nearThis, farThis);

    const scalar minDistThis =
        sqr(nearThis.x() - pt.x())
     +  sqr(nearThis.y() - pt.y())
     +  sqr(nearThis.z() - pt.z());
    const scalar maxDistThis =
        sqr(farThis.x() - pt.x())
     +  sqr(farThis.y() - pt.y())
     +  sqr(farThis.z() - pt.z());

    //
    // Distance point <-> other
    //

    point nearOther, farOther;

    // get nearest and furthest away vertex
    other.calcExtremities(pt, nearOther, farOther);

    const scalar minDistOther =
        sqr(nearOther.x() - pt.x())
     +  sqr(nearOther.y() - pt.y())
     +  sqr(nearOther.z() - pt.z());
    const scalar maxDistOther =
        sqr(farOther.x() - pt.x())
     +  sqr(farOther.y() - pt.y())
     +  sqr(farOther.z() - pt.z());

    //
    // Categorize
    //
    if (maxDistThis < minDistOther)
    {
        // All vertices of this are nearer to point than any vertex of other
        return -1;
    }
    else if (minDistThis > maxDistOther)
    {
        // All vertices of this are further from point than any vertex of other
        return 1;
    }
    else
    {
        // Mixed bag
        return 0;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

bool Foam::operator==(const treeBoundBox& a, const treeBoundBox& b)
{
    return operator==
    (
        static_cast<const boundBox&>(a),
        static_cast<const boundBox&>(b)
    );
}


bool Foam::operator!=(const treeBoundBox& a, const treeBoundBox& b)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const treeBoundBox& bb)
{
    return os << static_cast<const boundBox&>(bb);
}


Foam::Istream& Foam::operator>>(Istream& is, treeBoundBox& bb)
{
    return is >> static_cast<boundBox&>(bb);
}


// ************************************************************************* //
