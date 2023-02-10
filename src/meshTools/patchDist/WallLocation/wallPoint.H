/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::wallPoint

Description
    Holds information regarding nearest wall point. Used in wall distance
    calculation.

SourceFiles
    wallPointI.H

\*---------------------------------------------------------------------------*/

#ifndef wallPoint_H
#define wallPoint_H

#include "pointField.H"
#include "face.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class transformer;

// Forward declaration of friend functions and operators
class wallPoint;
Ostream& operator<<(Ostream&, const wallPoint&);
Istream& operator>>(Istream&, wallPoint&);

/*---------------------------------------------------------------------------*\
                          Class wallPoint Declaration
\*---------------------------------------------------------------------------*/

class wallPoint
{
    // Private Data

        //- Position of nearest wall center
        point origin_;

        //- Normal distance (squared) from cellcenter to origin
        scalar distSqr_;


protected:

    // Protected Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer pt. Return true if w2 is closer to point, false
        //  otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point&,
            const wallPoint& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    // Constructors

        //- Construct null
        inline wallPoint();

        //- Construct from origin, distance
        inline wallPoint(const point& origin, const scalar distSqr);

        //- Construct from face, distance
        inline wallPoint
        (
            const face& f,
            const pointField& ps,
            const point& centre,
            const scalar distSqr
        );


    // Member Functions

        // Access

            inline const point& origin() const;

            inline point& origin();

            inline scalar distSqr() const;

            inline scalar& distSqr();

            template<class TrackingData>
            inline scalar dist(TrackingData& td) const;


        // Needed by FaceCellWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data. Used for cyclics checking.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const wallPoint&,
                const scalar,
                TrackingData& td
            ) const;

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const transformer& transform,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal
            (
                const wallPoint&,
                TrackingData& td
            ) const;


    // Member Operators

        inline bool operator==(const wallPoint&) const;
        inline bool operator!=(const wallPoint&) const;


    // IOstream Operators

        inline friend Ostream& operator<<(Ostream&, const wallPoint&);
        inline friend Istream& operator>>(Istream&, wallPoint&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallPointI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
