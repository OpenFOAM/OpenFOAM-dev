/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
#include "Random.H"
#include "FixedList.H"
#include "UList.H"
#include "triPointRef.H"
#include "boundBox.H"
#include "barycentric.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

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

/*---------------------------------------------------------------------------*\
                         class tetrahedron Declaration
\*---------------------------------------------------------------------------*/

template<class Point, class PointRef>
class tetrahedron
{
    // Private Data

        PointRef a_, b_, c_, d_;


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
            inline triPointRef tri(const label facei) const;

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

            //- Calculate the point from the given barycentric coordinates.
            inline Point barycentricToPoint(const barycentric& bary) const;

            //- Calculate the barycentric coordinates from the given point
            inline barycentric pointToBarycentric(const point& pt) const;

            //- Calculate the barycentric coordinates from the given point.
            //  Returns the determinant.
            inline scalar pointToBarycentric
            (
                const point& pt,
                barycentric& bary
            ) const;

            //- Return nearest point to p on tetrahedron. Is p itself
            //  if inside.
            inline pointHit nearestPoint(const point& p) const;

            //- Return true if point is inside tetrahedron
            inline bool inside(const point& pt) const;

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

            //- Calculate the bounding box
            boundBox bounds() const;


    // IOstream Operators

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
    #include "tetrahedron.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
