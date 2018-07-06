/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::externalPointEdgePoint

Description
    Holds information regarding nearest wall point. Used in PointEdgeWave.
    (so not standard FaceCellWave)
    To be used in wall distance calculation.

SourceFiles
    externalPointEdgePointI.H
    externalPointEdgePoint.C

\*---------------------------------------------------------------------------*/

#ifndef externalPointEdgePoint_H
#define externalPointEdgePoint_H

#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;


// Forward declaration of friend functions and operators

class externalPointEdgePoint;

Istream& operator>>(Istream&, externalPointEdgePoint&);
Ostream& operator<<(Ostream&, const externalPointEdgePoint&);


/*---------------------------------------------------------------------------*\
                           Class externalPointEdgePoint Declaration
\*---------------------------------------------------------------------------*/

class externalPointEdgePoint
{
    // Private data

        //- Position of nearest wall center
        point origin_;

        //- Normal distance (squared) from point to origin
        scalar distSqr_;


    // Private Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer pt. Return true if w2 is closer to point,
        //  false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point&,
            const externalPointEdgePoint& w2,
            const scalar tol,
            TrackingData& td
        );

        //- Combine current with w2. Update distSqr, origin if w2 has smaller
        //  quantities and returns true.
        template<class TrackingData>
        inline bool update
        (
            const externalPointEdgePoint& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    //- Class used to pass data into container
    class trackingData
    {
    public:
        const pointField& points_;

        trackingData(const pointField& points)
        :
            points_(points)
        {}
    };



    // Constructors

        //- Construct null
        inline externalPointEdgePoint();

        //- Construct from origin, distance
        inline externalPointEdgePoint(const point&, const scalar);

        //- Construct as copy
        inline externalPointEdgePoint(const externalPointEdgePoint&);


    // Member Functions

        // Access

            inline const point& origin() const;

            inline scalar distSqr() const;


        // Needed by PointEdgeWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data. Used for cyclics checking.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const externalPointEdgePoint&,
                const scalar tol,
                TrackingData& td
            ) const;

            //- Convert origin to relative vector to leaving point
            //  (= point coordinate)
            template<class TrackingData>
            inline void leaveDomain
            (
                const polyPatch& patch,
                const label patchPointi,
                const point& pos,
                TrackingData& td
            );

            //- Convert relative origin to absolute by adding entering point
            template<class TrackingData>
            inline void enterDomain
            (
                const polyPatch& patch,
                const label patchPointi,
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
                const label pointi,
                const label edgeI,
                const externalPointEdgePoint& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  Merge new and old info.
            template<class TrackingData>
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointi,
                const externalPointEdgePoint& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  No information about current position whatsoever.
            template<class TrackingData>
            inline bool updatePoint
            (
                const externalPointEdgePoint& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of point on edge.
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const label edgeI,
                const label pointi,
                const externalPointEdgePoint& pointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Equivalent to operator== with TrackingData
            template<class TrackingData>
            inline bool equal
            (
                const externalPointEdgePoint&,
                TrackingData& td
            ) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const externalPointEdgePoint&) const;
        inline bool operator!=(const externalPointEdgePoint&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const externalPointEdgePoint&);
        friend Istream& operator>>(Istream&, externalPointEdgePoint&);
};


//- Data associated with externalPointEdgePoint type are contiguous
template<>
inline bool contiguous<externalPointEdgePoint>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "externalPointEdgePointI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
