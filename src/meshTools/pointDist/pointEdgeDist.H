/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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
    Foam::pointEdgeDist

Description
    Holds information regarding nearest wall point. Used in PointEdgeWave.
    (so not standard FaceCellWave)
    To be used in wall distance calculation.

SourceFiles
    pointEdgeDistI.H
    pointEdgeDist.C

\*---------------------------------------------------------------------------*/

#ifndef pointEdgeDist_H
#define pointEdgeDist_H

#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;
class transformer;


// Forward declaration of friend functions and operators

class pointEdgeDist;

Istream& operator>>(Istream&, pointEdgeDist&);
Ostream& operator<<(Ostream&, const pointEdgeDist&);


/*---------------------------------------------------------------------------*\
                   Class pointEdgeDist Declaration
\*---------------------------------------------------------------------------*/

class pointEdgeDist
{
    // Private Data

        //- Position of nearest wall center
        point origin_;

        //- Normal distance (squared) from point to origin
        scalar distSqr_;


public:

    //- Class used to pass data into container
    class data
    {
    public:
        const pointField& points;

        data(const pointField& p)
        :
            points(p)
        {}
    };


private:

    // Private Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer pt. Return true if w2 is closer to point,
        //  false otherwise.
        inline bool update
        (
            const point&,
            const pointEdgeDist& w2,
            const scalar tol,
            data& td
        );


public:

    // Constructors

        //- Construct null
        inline pointEdgeDist();

        //- Construct from origin, distance
        inline pointEdgeDist(const point&, const scalar);


    // Member Functions

        // Access

            inline const point& origin() const;

            inline scalar distSqr() const;


        // Needed by PointEdgeWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            inline bool valid(data& td) const;

            //- Transform across an interface
            inline void transform
            (
                const polyPatch& patch,
                const label patchFacei,
                const transformer& transform,
                data& td
            );

            //- Influence of edge on point
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointi,
                const label edgei,
                const pointEdgeDist& edgeinfo,
                const scalar tol,
                data& td
            );

            //- Influence of different value on same point.
            //  Merge new and old info.
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointi,
                const pointEdgeDist& newPointInfo,
                const scalar tol,
                data& td
            );

            //- Influence of point on edge.
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const label edgei,
                const label pointi,
                const pointEdgeDist& pointInfo,
                const scalar tol,
                data& td
            );

            //- Equivalent to operator== with data
            inline bool equal
            (
                const pointEdgeDist&,
                data& td
            ) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const pointEdgeDist&) const;
        inline bool operator!=(const pointEdgeDist&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const pointEdgeDist&);
        friend Istream& operator>>(Istream&, pointEdgeDist&);
};


//- Data associated with pointEdgeDist type are contiguous
template<>
inline bool contiguous<pointEdgeDist>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointEdgeDistI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
