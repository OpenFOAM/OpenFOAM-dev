/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::boundSphere

Description
    The smallest sphere enclosing a given set of points

SourceFiles
    boundSphereI.H
    boundSphereTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef boundSphere_H
#define boundSphere_H

#include "point.H"
#include "DLPermutation.H"
#include "nil.H"
#include "RemoteData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class boundSphere;

bool operator==(const boundSphere&, const boundSphere&);
bool operator!=(const boundSphere&, const boundSphere&);

Istream& operator>>(Istream&, boundSphere&);
Ostream& operator<<(Ostream&, const boundSphere&);

/*---------------------------------------------------------------------------*\
                          Class boundSphere Declaration
\*---------------------------------------------------------------------------*/

class boundSphere
{
    // Private Data

        //- The centre
        point c_;

        //- The radius squared
        scalar rSqr_;


    // Private Type Definitions

        //- For a given List type, return the corresponding UList-type
        template<class PointList>
        using PointUList =
            typename std::conditional
            <
                std::is_base_of<UList<point>, PointList>::value,
                UList<point>,
                typename std::conditional
                <
                    std::is_base_of<UIndirectList<point>, PointList>::value,
                    UIndirectList<point>,
                    nil
                >::type
            >::type;


    // Private Classes

        //- A structure to manage a fixed set of local points and a growing
        //  list of added remote points and provides indexing into these sets
        //  as if they were in a single contiguous list. It is used when
        //  constructing global bound spheres. A global bound sphere is
        //  calculated by first constructing the local bound sphere and then
        //  synchronising the limiting set of outermost points. These points
        //  are then added to the list and the process is repeated until
        //  no more points are added and all processes arrive at the same
        //  set of boundary points.
        template<class PointList>
        class LocalAndRemotePoints
        {
            const PointList& locals_;

            DynamicList<RemoteData<point>> remotes_;

            HashTable<label, remote, Hash<remote>> remoteTable_;

        public:

            inline LocalAndRemotePoints(const PointList& locals)
            :
                locals_(locals)
            {}

            inline label nRemotes() const
            {
                return remotes_.size();
            }

            inline const point& operator[](const label i) const
            {
                return
                    i < locals_.size()
                  ? locals_[i]
                  : remotes_[i - locals_.size()].data;
            }

            inline RemoteData<point> operator()(const label i) const
            {
                return
                    i < locals_.size()
                  ? RemoteData<point>(Pstream::myProcNo(), i, locals_[i])
                  : remotes_[i - locals_.size()];
            }

            inline bool found(const remote& r) const
            {
                if (r.proci == Pstream::myProcNo()) return true;
                return remoteTable_.found(r);
            }

            inline label insert(const RemoteData<point>& rp)
            {
                if (rp.proci == Pstream::myProcNo()) return rp.elementi;
                if (remoteTable_.found(rp)) return remoteTable_[rp];
                return insertNoCheck(rp);
            }

            inline label insertNoCheck(const RemoteData<point>& rp)
            {
                remotes_.append(rp);
                remoteTable_.insert(rp, remotes_.size() - 1);
                return locals_.size() + remotes_.size();
            }
        };


    // Private Constructor

        //- Construct from components
        inline boundSphere(const point& c, const scalar rSqr);


    // Private Static Member Functions

        //- Cast a List to the corresponding UList-type
        template<class PointList>
        static inline PointUList<PointList> pointUList(const PointList&);

        //- Sort boundary point indices
        static inline void sortBoundaryPis(FixedList<label, 4>&);

        //- Sort boundary remote point indices
        template<class Remote>
        static inline void sortBoundaryPis(FixedList<Remote, 4>&);

        //- ...
        struct implementation
        {

        //- Return the sphere intersecting the given set of up to four points
        template<class PointList>
        static boundSphere intersect
        (
            const PointList& ps,
            const FixedList<label, 4>& pis,
            const label nPs
        );

        //- Return the sphere bounding the given set of up to four points
        template<class PointList>
        static boundSphere trivial
        (
            const PointList& ps,
            const FixedList<label, 4>& pis,
            const label nPs,
            FixedList<label, 4>& boundaryPis
        );

        //- Return the sphere bounding the given set of points using an
        //  inefficient brute-force method
        template<class PointList>
        static boundSphere bruteForce
        (
            const PointList& ps,
            FixedList<label, 4>& boundaryPis
        );

        //- Inner/recursive implementation of Welzl's algorithm
        template<class PointList, class DLPermutation>
        static boundSphere Welzl
        (
            const PointList& ps,
            DLPermutation& pis,
            const typename DLPermutation::const_iterator& pisEnd,
            FixedList<label, 4>& boundaryPis,
            const label nBoundaryPs
        );

        //- Return the sphere bounding the given set of points using Welzl's
        //  algorithm. Constructs the local sphere containing the points on
        //  this processor only.
        template<class PointList>
        static boundSphere local
        (
            const PointList& ps,
            randomGenerator& rndGen,
            FixedList<label, 4>& boundaryPis
        );

        //- Return the sphere bounding the given set of points using Welzl's
        //  algorithm. Constructs the global sphere containing the points on
        //  all processors.
        template<class PointList>
        static boundSphere global
        (
            const PointList& ps,
            randomGenerator& rndGen,
            FixedList<RemoteData<point>, 4>& boundaryPis,
            const bool strict
        );

        //- ...
        };


public:


    // Constructors

        //- Construct a null/invalid sphere
        inline boundSphere();

        //- Construct given a point and a radius
        static inline boundSphere cr(const point& c, const scalar r);

        //- Construct given a point and a radius
        static inline boundSphere cr(const Tuple2<point, scalar>& cr);

        //- Construct given a point and a radius squared
        static inline boundSphere crSqr(const point& c, const scalar rSqr);

        //- Construct given a point and a radius squared
        static inline boundSphere crSqr(const Tuple2<point, scalar>& crSqr);


    // Member Functions

        //- Return the centre
        inline const point& c() const;

        //- Return the radius squared
        inline scalar rSqr() const;

        //- Return the radius
        inline scalar r() const;

        //- Return whether this sphere is valid
        inline bool valid() const;

        //- Return whether this sphere contains a given point
        inline bool contains(const point& p) const;

        //- Expand the sphere by the given factor
        inline void inflate(const scalar s);


    // Static Member Functions

        //- Return whether two spheres overlap
        static inline bool overlap(const boundSphere& a, const boundSphere& b);

        //- Return the sphere intersecting the given set of up to four points
        template<class PointList>
        static inline boundSphere intersect
        (
            const PointList& ps,
            const FixedList<label, 4>& pis,
            const label nPs
        );

        //- Return the sphere bounding the given set of up to four points
        template<class PointList>
        static inline boundSphere trivial
        (
            const PointList& ps,
            const FixedList<label, 4>& pis,
            const label nPs
        );

        //- Return the sphere bounding the given set of up to four points
        template<class PointList>
        static inline boundSphere trivial
        (
            const PointList& ps,
            const FixedList<label, 4>& pis,
            const label nPs,
            FixedList<label, 4>& boundaryPis
        );

        //- Return the sphere bounding the given set of points using an
        //  inefficient brute-force method
        template<class PointList>
        static inline boundSphere bruteForce
        (
            const PointList& ps
        );

        //- Return the sphere bounding the given set of points using an
        //  inefficient brute-force method
        template<class PointList>
        static inline boundSphere bruteForce
        (
            const PointList& ps,
            FixedList<label, 4>& boundaryPis
        );

        //- Return the sphere bounding the given set of points using Welzl's
        //  algorithm. Constructs the local sphere containing the points on
        //  this processor only.
        template<class PointList>
        static inline boundSphere local
        (
            const PointList& ps
        );

        //- Return the sphere bounding the given set of points using Welzl's
        //  algorithm. Constructs the local sphere containing the points on
        //  this processor only.
        template<class PointList>
        static inline boundSphere local
        (
            const PointList& ps,
            randomGenerator& rndGen
        );

        //- Return the sphere bounding the given set of points using Welzl's
        //  algorithm. Constructs the local sphere containing the points on
        //  this processor only.
        template<class PointList>
        static inline boundSphere local
        (
            const PointList& ps,
            FixedList<label, 4>& boundaryPis
        );

        //- Return the sphere bounding the given set of points using Welzl's
        //  algorithm. Constructs the local sphere containing the points on
        //  this processor only.
        template<class PointList>
        static inline boundSphere local
        (
            const PointList& ps,
            randomGenerator& rndGen,
            FixedList<label, 4>& boundaryPis
        );

        //- Return the sphere bounding the given set of points using Welzl's
        //  algorithm. Constructs the global sphere containing the points on
        //  all processors.
        template<class PointList>
        static inline boundSphere global
        (
            const PointList& ps,
            const bool strict = false
        );

        //- Return the sphere bounding the given set of points using Welzl's
        //  algorithm. Constructs the global sphere containing the points on
        //  all processors.
        template<class PointList>
        static inline boundSphere global
        (
            const PointList& ps,
            randomGenerator& rndGen,
            const bool strict = false
        );

        //- Return the sphere bounding the given set of points using Welzl's
        //  algorithm. Constructs the global sphere containing the points on
        //  all processors.
        template<class PointList>
        static inline boundSphere global
        (
            const PointList& ps,
            FixedList<RemoteData<point>, 4>& boundaryPis,
            const bool strict = false
        );

        //- Return the sphere bounding the given set of points using Welzl's
        //  algorithm. Constructs the global sphere containing the points on
        //  all processors.
        template<class PointList>
        static inline boundSphere global
        (
            const PointList& ps,
            randomGenerator& rndGen,
            FixedList<RemoteData<point>, 4>& boundaryPis,
            const bool strict = false
        );


    // Friend Operators

        //- Equality comparison
        inline friend bool operator==(const boundSphere&, const boundSphere&);

        //- Inequality comparison
        inline friend bool operator!=(const boundSphere&, const boundSphere&);


    // IOstream operators

        //- Read form stream
        inline friend Istream& operator>>(Istream&, boundSphere&);

        //- Write to stream
        inline friend Ostream& operator<<(Ostream&, const boundSphere&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "boundSphereI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
