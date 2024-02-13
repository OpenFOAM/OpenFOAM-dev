/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::pistonPointEdgeData

Description
    Holds information regarding nearest wall point. Used in PointEdgeWave.
    (so not standard FaceCellWave)
    To be used in wall distance calculation.

SourceFiles
    pistonPointEdgeDataI.H
    pistonPointEdgeData.C

\*---------------------------------------------------------------------------*/

#ifndef pistonPointEdgeData_H
#define pistonPointEdgeData_H

#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;
class transformer;


// Forward declaration of friend functions and operators

class pistonPointEdgeData;

Istream& operator>>(Istream&, pistonPointEdgeData&);
Ostream& operator<<(Ostream&, const pistonPointEdgeData&);


/*---------------------------------------------------------------------------*\
                   Class pistonPointEdgeData Declaration
\*---------------------------------------------------------------------------*/

class pistonPointEdgeData
{
    // Private Data

        //- Has this point or edge been visited
        bool visited_;

        //- Is this point in the bowl
        bool inBowl_;


public:

    //- Class used to pass data into container
    class trackingData
    {
    public:
        const pointField& points_;
        const vector axis_;
        const scalar maxZ_;

        trackingData
        (
            const pointField& points,
            const vector& axis,
            const scalar maxZ
        )
        :
            points_(points),
            axis_(axis),
            maxZ_(maxZ)
        {}
    };


private:

    // Private Member Functions

        //- Update visited and inBowl
        inline bool update(const point&, trackingData& td);


public:

    // Constructors

        //- Construct null
        inline pistonPointEdgeData();

        //- Construct from inBowl
        inline pistonPointEdgeData(const bool inBowl);


    // Member Functions

        // Access

            inline bool inBowl() const;


        // Needed by PointEdgeWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            inline bool valid(trackingData& td) const;

            //- Transform across an interface
            inline void transform
            (
                const polyPatch& patch,
                const label patchFacei,
                const transformer& transform,
                trackingData& td
            );

            //- Influence of edge on point
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointi,
                const label edgei,
                const pistonPointEdgeData& edgeinfo,
                const scalar tol,
                trackingData& td
            );

            //- Influence of different value on same point.
            //  Merge new and old info.
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointi,
                const pistonPointEdgeData& newPointInfo,
                const scalar tol,
                trackingData& td
            );

            //- Influence of point on edge.
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const label edgei,
                const label pointi,
                const pistonPointEdgeData& pointInfo,
                const scalar tol,
                trackingData& td
            );

            //- Equivalent to operator== with trackingData
            inline bool equal
            (
                const pistonPointEdgeData&,
                trackingData& td
            ) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const pistonPointEdgeData&) const;
        inline bool operator!=(const pistonPointEdgeData&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const pistonPointEdgeData&);
        friend Istream& operator>>(Istream&, pistonPointEdgeData&);
};


//- Data associated with pistonPointEdgeData type are contiguous
template<>
inline bool contiguous<pistonPointEdgeData>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pistonPointEdgeDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
