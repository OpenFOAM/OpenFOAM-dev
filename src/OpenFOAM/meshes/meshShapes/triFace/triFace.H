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

Class
    Foam::triFace

Description
    A triangular face using a FixedList of labels corresponding to mesh
    vertices.

See also
    Foam::face, Foam::triangle

SourceFiles
    triFaceI.H
    triFaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef triFace_H
#define triFace_H

#include "FixedList.H"
#include "edgeList.H"
#include "pointHit.H"
#include "intersection.H"
#include "pointField.H"
#include "triPointRef.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class face;
class triFace;

inline bool operator==(const triFace&, const triFace&);
inline bool operator!=(const triFace&, const triFace&);


/*---------------------------------------------------------------------------*\
                          Class triFace Declaration
\*---------------------------------------------------------------------------*/

class triFace
:
    public FixedList<label, 3>
{

public:

    // Constructors

        //- Construct null
        inline triFace();

        //- Construct from three point labels
        inline triFace
        (
            const label a,
            const label b,
            const label c
        );

        //- Construct from a list of labels
        explicit inline triFace(const labelUList&);

        //- Construct from Istream
        inline triFace(Istream&);


    // Member Functions

        //- Collapse face by removing duplicate point labels
        //  return the collapsed size, set collapsed point labels to -1
        inline label collapse();

        //- Flip the face in-place.
        //  The starting points of the original and flipped face are identical.
        inline void flip();

        //- Return the points corresponding to this face
        inline pointField points(const pointField&) const;

        //- Return triangle as a face
        inline face triFaceFace() const;

        //- Return the triangle
        inline triPointRef tri(const pointField&) const;

        //- Return centre (centroid)
        inline point centre(const pointField&) const;

        //- Calculate average value at centroid of face
        template<class Type>
        Type average(const pointField&, const Field<Type>&) const;

        //- Return vector area
        inline vector area(const pointField&) const;

        //- Return scalar magnitude
        inline scalar mag(const pointField&) const;

        //- Return unit normal
        inline vector normal(const pointField&) const;

        //- Number of triangles after splitting
        inline label nTriangles() const;

        //- Return face with reverse direction
        //  The starting points of the original and reverse face are identical.
        inline triFace reverseFace() const;

        //- Return swept-volume
        inline scalar sweptVol
        (
            const pointField& oldPoints,
            const pointField& newPoints
        ) const;

        //- Return the inertia tensor, with optional reference
        //  point and density specification
        inline tensor inertia
        (
            const pointField&,
            const point& refPt = vector::zero,
            scalar density = 1.0
        ) const;

        //- Return point intersection with a ray starting at p,
        //  with direction q.
        inline pointHit ray
        (
            const point& p,
            const vector& q,
            const pointField& points,
            const intersection::algorithm = intersection::FULL_RAY,
            const intersection::direction dir = intersection::VECTOR
        ) const;

        //- Fast intersection with a ray.
        inline pointHit intersection
        (
            const point& p,
            const vector& q,
            const pointField& points,
            const intersection::algorithm alg,
            const scalar tol = 0.0
        ) const;

        inline pointHit intersection
        (
            const point& p,
            const vector& q,
            const point& ctr,
            const pointField& points,
            const intersection::algorithm alg,
            const scalar tol = 0.0
        ) const;

        //- Return nearest point to face
        inline pointHit nearestPoint
        (
            const point& p,
            const pointField& points
        ) const;


        //- Return nearest point to face and classify it:
        //  + near point (nearType=POINT, nearLabel=0, 1, 2)
        //  + near edge (nearType=EDGE, nearLabel=0, 1, 2)
        //    Note: edges are counted from starting vertex so
        //    e.g. edge n is from f[n] to f[0], where the face has n + 1
        //    points
        inline pointHit nearestPointClassify
        (
            const point& p,
            const pointField& points,
            label& nearType,
            label& nearLabel
        ) const;

        //- Return number of edges
        inline label nEdges() const;

        //- Return edges in face point ordering,
        //  i.e. edges()[0] is edge between [0] and [1]
        inline edgeList edges() const;

        //- Return n-th face edge
        inline edge faceEdge(const label n) const;

        //- Return the edge direction on the face
        //  Returns:
        //  - +1: forward (counter-clockwise) on the face
        //  - -1: reverse (clockwise) on the face
        //  -  0: edge not found on the face
        inline int edgeDirection(const edge&) const;

        //- Compare triFaces
        //  Returns:
        //  -  0: different
        //  - +1: identical
        //  - -1: same face, but different orientation
        static inline int compare(const triFace&, const triFace&);

    // Friend Operators

        inline friend bool operator==(const triFace&, const triFace&);
        inline friend bool operator!=(const triFace&, const triFace&);
};


//- Hash specialization for hashing triFace - a commutative hash value.
//  Hash incrementally.
template<>
inline unsigned Hash<triFace>::operator()(const triFace& t, unsigned seed) const
{
    // Fortunately we don't need this very often
    const uLabel t0(t[0]);
    const uLabel t1(t[1]);
    const uLabel t2(t[2]);

    const uLabel val = (t0*t1*t2 + t0+t1+t2);

    return Hash<uLabel>()(val, seed);
}


//- Hash specialization for hashing triFace - a commutative hash value.
template<>
inline unsigned Hash<triFace>::operator()(const triFace& t) const
{
    return Hash<triFace>::operator()(t, 0);
}


template<>
inline bool contiguous<triFace>()  {return true;}


//- Hash specialization to offset faces in ListListOps::combineOffset
template<>
class offsetOp<triFace>
{

public:

    inline triFace operator()
    (
        const triFace& x,
        const label offset
    ) const
    {
        triFace result(x);

        forAll(x, xI)
        {
            result[xI] = x[xI] + offset;
        }
        return result;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triFaceI.H"

#ifdef NoRepository
    #include "triFaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
