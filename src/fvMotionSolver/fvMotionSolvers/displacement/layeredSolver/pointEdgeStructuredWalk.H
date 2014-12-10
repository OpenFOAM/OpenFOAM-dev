/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::pointEdgeStructuredWalk

Description
    Determines length of string of edges walked to point.

SourceFiles
    pointEdgeStructuredWalkI.H
    pointEdgeStructuredWalk.C

\*---------------------------------------------------------------------------*/

#ifndef pointEdgeStructuredWalk_H
#define pointEdgeStructuredWalk_H

#include "point.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class pointEdgeStructuredWalk Declaration
\*---------------------------------------------------------------------------*/

class pointEdgeStructuredWalk
{
    // Private data

        //- Starting location
        point point0_;

        //- Previous point
        point previousPoint_;

        //- Sum of distance
        scalar dist_;

        //- Passive data
        vector data_;

    // Private Member Functions

        //- Evaluate distance to point.
        template<class TrackingData>
        inline bool update
        (
            const pointEdgeStructuredWalk& w2,
            const scalar tol,
            TrackingData& td
        );

public:

    // Constructors

        //- Construct null
        inline pointEdgeStructuredWalk();

        //- Construct from components
        inline pointEdgeStructuredWalk
        (
            const point&,
            const point&,
            const scalar,
            const vector&
        );


    // Member Functions

        // Access

            inline bool inZone() const;

            inline scalar dist() const;

            inline const vector& data() const;


        // Needed by meshWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data. Used for cyclics checking.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const pointEdgeStructuredWalk&,
                const scalar tol,
                TrackingData& td
            ) const;

            //- Convert origin to relative vector to leaving point
            //  (= point coordinate)
            template<class TrackingData>
            inline void leaveDomain
            (
                const polyPatch& patch,
                const label patchPointI,
                const point& pos,
                TrackingData& td
            );

            //- Convert relative origin to absolute by adding entering point
            template<class TrackingData>
            inline void enterDomain
            (
                const polyPatch& patch,
                const label patchPointI,
                const point& pos,
                TrackingData& td
            );

            //- Apply rotation matrix to origin
            template<class TrackingData>
            inline void transform
            (
                const tensor& rotTensor,
                TrackingData& td
            );

            //- Influence of edge on point
            template<class TrackingData>
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointI,
                const label edgeI,
                const pointEdgeStructuredWalk& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  Merge new and old info.
            template<class TrackingData>
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointI,
                const pointEdgeStructuredWalk& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  No information about current position whatsoever.
            template<class TrackingData>
            inline bool updatePoint
            (
                const pointEdgeStructuredWalk& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of point on edge.
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const label edgeI,
                const label pointI,
                const pointEdgeStructuredWalk& pointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal(const pointEdgeStructuredWalk&, TrackingData&)
            const;


    // Member Operators

        //Note: Used to determine whether to call update.
        inline bool operator==(const pointEdgeStructuredWalk&) const;
        inline bool operator!=(const pointEdgeStructuredWalk&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const pointEdgeStructuredWalk&);
        friend Istream& operator>>(Istream&, pointEdgeStructuredWalk&);
};


//- Data associated with pointEdgeStructuredWalk type are contiguous
template<>
inline bool contiguous<pointEdgeStructuredWalk>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointEdgeStructuredWalkI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
