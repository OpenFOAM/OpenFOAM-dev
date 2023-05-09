/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::pointEdgeLayerInfo

Description
    Class to be used with PointEdgeWave which enumerates layers of points

SourceFiles
    pointEdgeLayerInfoI.H

\*---------------------------------------------------------------------------*/

#ifndef pointEdgeLayerInfo_H
#define pointEdgeLayerInfo_H

#include "contiguous.H"
#include "label.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;
class transformer;

// Forward declaration of friend functions and operators
class pointEdgeLayerInfo;
inline Ostream& operator<<(Ostream&, const pointEdgeLayerInfo&);
inline Istream& operator>>(Istream&, pointEdgeLayerInfo&);

/*---------------------------------------------------------------------------*\
                          Class pointEdgeLayerInfo Declaration
\*---------------------------------------------------------------------------*/

class pointEdgeLayerInfo
{
    // Private Data

        //- Layer index
        label layer_;


public:

    // Constructors

        //- Construct null
        inline pointEdgeLayerInfo();

        //- Construct given the point layer index
        inline pointEdgeLayerInfo(const label pointLayer);


    // Member Functions

        // Access

            //- Return the point layer index
            inline label pointLayer() const;

            //- Return the edge layer index
            inline label edgeLayer() const;


        // Needed by PointEdgeWave

            //- Check whether info has been changed at all or
            //  still contains original (invalid) value
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data. Used for cyclics checking.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const pointEdgeLayerInfo&,
                const scalar tol,
                TrackingData& td
            ) const;

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const polyPatch& patch,
                const label patchPointi,
                const transformer& transform,
                TrackingData& td
            );

            //- Influence of edge on point
            template<class TrackingData>
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointi,
                const label edgei,
                const pointEdgeLayerInfo& edgeInfo,
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
                const pointEdgeLayerInfo& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  No old information.
            template<class TrackingData>
            inline bool updatePoint
            (
                const pointEdgeLayerInfo& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of point on edge
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const label edgei,
                const label pointi,
                const pointEdgeLayerInfo& pointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal
            (
                const pointEdgeLayerInfo&,
                TrackingData& td
            ) const;


    // Member Operators

        inline bool operator==(const pointEdgeLayerInfo&) const;
        inline bool operator!=(const pointEdgeLayerInfo&) const;


    // IOstream Operators

        inline friend Ostream& operator<<(Ostream&, const pointEdgeLayerInfo&);
        inline friend Istream& operator>>(Istream&, pointEdgeLayerInfo&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
inline bool contiguous<pointEdgeLayerInfo>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointEdgeLayerInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
