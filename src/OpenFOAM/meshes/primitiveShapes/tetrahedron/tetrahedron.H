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

Class
    Foam::tetrahedron

Description
    A tetrahedron primitive.

    Ordering of edges needs to be the same for a tetrahedron
    class, a tetrahedron cell shape model and a tetCell.

SourceFiles
    tetrahedronI.H
    tetrahedron.C

\*---------------------------------------------------------------------------*/

#ifndef tetrahedron_H
#define tetrahedron_H

#include "point.H"
#include "primitiveFieldsFwd.H"
#include "pointHit.H"
#include "cachedRandom.H"
#include "Random.H"
#include "FixedList.H"
#include "UList.H"
#include "triPointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;
class tetPoints;
class plane;

// Forward declaration of friend functions and operators

template<class Point, class PointRef> class tetrahedron;

template<class Point, class PointRef>
inline Istream& operator>>
(
    Istream&,
    tetrahedron<Point, PointRef>&
);

template<class Point, class PointRef>
inline Ostream& operator<<
(
    Ostream&,
    const tetrahedron<Point, PointRef>&
);

typedef tetrahedron<point, const point&> tetPointRef;

/*---------------------------------------------------------------------------*\
                           class tetrahedron Declaration
\*---------------------------------------------------------------------------*/

template<class Point, class PointRef>
class tetrahedron
{
public:

    // Public typedefs

        //- Storage type for tets originating from intersecting tets.
        //  (can possibly be smaller than 200)
        typedef FixedList<tetPoints, 200> tetIntersectionList;


        // Classes for use in sliceWithPlane. What to do with decomposition
        // of tet.

            //- Dummy
            class dummyOp
            {
            public:
                inline void operator()(const tetPoints&);
            };

            //- Sum resulting volumes
            class sumVolOp
            {
            public:
                scalar vol_;

                inline sumVolOp();

                inline void operator()(const tetPoints&);
            };

            //- Store resulting tets
            class storeOp
            {
                tetIntersectionList& tets_;
                label& nTets_;

            public:
                inline storeOp(tetIntersectionList&, label&);

                inline void operator()(const tetPoints&);
            };

private:

    // Private data

        PointRef a_, b_, c_, d_;

        inline static point planeIntersection
        (
            const FixedList<scalar, 4>&,
            const tetPoints&,
            const label,
            const label
        );

        template<class TetOp>
        inline static void decomposePrism
        (
            const FixedList<point, 6>& points,
            TetOp& op
        );

        template<class AboveTetOp, class BelowTetOp>
        inline static void tetSliceWithPlane
        (
            const plane& pl,
            const tetPoints& tet,
            AboveTetOp& aboveOp,
            BelowTetOp& belowOp
        );


public:

    // Member constants

        enum
        {
            nVertices = 4,  // Number of vertices in tetrahedron
            nEdges = 6      // Number of edges in tetrahedron
        };


    // Constructors

        //- Construct from points
        inline tetrahedron
        (
            const Point& a,
            const Point& b,
            const Point& c,
            const Point& d
        );

        //- Construct from four points in the list of points
        inline tetrahedron
        (
            const UList<Point>&,
            const FixedList<label, 4>& indices
        );

        //- Construct from Istream
        inline tetrahedron(Istream&);


    // Member Functions

        // Access

            //- Return vertices
            inline const Point& a() const;

            inline const Point& b() const;

            inline const Point& c() const;

            inline const Point& d() const;

            //- Return i-th face
            inline triPointRef tri(const label faceI) const;

        // Properties

            //- Return face normal
            inline vector Sa() const;

            inline vector Sb() const;

            inline vector Sc() const;

            inline vector Sd() const;

            //- Return centre (centroid)
            inline Point centre() const;

            //- Return volume
            inline scalar mag() const;

            //- Return circum-centre
            inline Point circumCentre() const;

            //- Return circum-radius
            inline scalar circumRadius() const;

            //- Return quality: Ratio of tetrahedron and circum-sphere
            //  volume, scaled so that a regular tetrahedron has a
            //  quality of 1
            inline scalar quality() const;

            //- Return a random point in the tetrahedron from a
            //  uniform distribution
            inline Point randomPoint(Random& rndGen) const;

            //- Return a random point in the tetrahedron from a
            //  uniform distribution
            inline Point randomPoint(cachedRandom& rndGen) const;

            //- Calculate the barycentric coordinates of the given
            //  point, in the same order as a, b, c, d.  Returns the
            //  determinant of the solution.
            inline scalar barycentric
            (
                const point& pt,
                List<scalar>& bary
            ) const;

            //- Return nearest point to p on tetrahedron. Is p itself
            //  if inside.
            inline pointHit nearestPoint(const point& p) const;

            //- Return true if point is inside tetrahedron
            inline bool inside(const point& pt) const;

            //- Decompose tet into tets above and below plane
            template<class AboveTetOp, class BelowTetOp>
            inline void sliceWithPlane
            (
                const plane& pl,
                AboveTetOp& aboveOp,
                BelowTetOp& belowOp
            ) const;

            //- Decompose tet into tets inside and outside other tet
            inline void tetOverlap
            (
                const tetrahedron<Point, PointRef>& tetB,
                tetIntersectionList& insideTets,
                label& nInside,
                tetIntersectionList& outsideTets,
                label& nOutside
            ) const;


            //- Return (min)containment sphere, i.e. the smallest sphere with
            //  all points inside. Returns pointHit with:
            //  - hit         : if sphere is equal to circumsphere
            //                  (biggest sphere)
            //  - point       : centre of sphere
            //  - distance    : radius of sphere
            //  - eligiblemiss: false
            // Tol (small compared to 1, e.g. 1e-9) is used to determine
            // whether point is inside: mag(pt - ctr) < (1+tol)*radius.
            pointHit containmentSphere(const scalar tol) const;

            //- Fill buffer with shape function products
            void gradNiSquared(scalarField& buffer) const;

            void gradNiDotGradNj(scalarField& buffer) const;

            void gradNiGradNi(tensorField& buffer) const;

            void gradNiGradNj(tensorField& buffer) const;


    // IOstream operators

        friend Istream& operator>> <Point, PointRef>
        (
            Istream&,
            tetrahedron&
        );

        friend Ostream& operator<< <Point, PointRef>
        (
            Ostream&,
            const tetrahedron&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tetrahedronI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "tetrahedron.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
