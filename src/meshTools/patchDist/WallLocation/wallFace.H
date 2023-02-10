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
    Foam::wallFace

Description
    Holds information regarding nearest wall point. Used in wall distance
    calculation.

SourceFiles
    wallFaceI.H

\*---------------------------------------------------------------------------*/

#ifndef wallFace_H
#define wallFace_H

#include "pointField.H"
#include "face.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class transformer;

// Forward declaration of friend functions and operators
class wallFace;
Ostream& operator<<(Ostream&, const wallFace&);
Istream& operator>>(Istream&, wallFace&);

/*---------------------------------------------------------------------------*\
                          Class wallFace Declaration
\*---------------------------------------------------------------------------*/

class wallFace
{
    // Private Data

        //- Face points
        pointField points_;

        //- Normal distance (squared) from cell center to face
        scalar distSqr_;


protected:

    // Protected Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer pt. Return true if w2 is closer to point, false
        //  otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point& pt,
            const wallFace& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    // Constructors

        //- Construct null
        inline wallFace();

        //- Construct from face, distance
        inline wallFace
        (
            const face& f,
            const pointField& points,
            const scalar distSqr
        );

        //- Construct from face, distance
        inline wallFace
        (
            const face& f,
            const pointField& points,
            const point& centre,
            const scalar distSqr
        );


    // Member Functions

        // Access

            inline const pointField& points() const;

            inline pointField& points();

            inline scalar distSqr() const;

            inline scalar& distSqr();

            template<class TrackingData>
            inline scalar dist(TrackingData& td) const;


        // Needed by FaceCellWave

            //- Check whether the wallFace has been changed at all or still
            //  contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data. Used for checking
            //  consistency across cyclics.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const wallFace&,
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

            //- Test equality
            template<class TrackingData>
            inline bool equal
            (
                const wallFace&,
                TrackingData& td
            ) const;


    // Member Operators

        inline bool operator==(const wallFace&) const;
        inline bool operator!=(const wallFace&) const;


    // IOstream Operators

        inline friend Ostream& operator<<(Ostream&, const wallFace&);
        inline friend Istream& operator>>(Istream&, wallFace&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallFaceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
