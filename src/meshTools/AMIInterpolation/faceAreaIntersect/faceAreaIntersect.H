/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::faceAreaIntersect

Description
    Face intersection class
    - calculates intersection area by sub-dividing face into triangles
      and cutting

SourceFiles
    faceAreaIntersect.C

\*---------------------------------------------------------------------------*/

#ifndef faceAreaIntersect_H
#define faceAreaIntersect_H

#include "pointField.H"
#include "FixedList.H"
#include "plane.H"
#include "face.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class faceAreaIntersect Declaration
\*---------------------------------------------------------------------------*/

class faceAreaIntersect
{
public:

    typedef FixedList<point, 3> triPoints;

    enum triangulationMode
    {
        tmFan,
        tmMesh
    };

    static const NamedEnum<triangulationMode, 2> triangulationModeNames_;


private:

    // Private data

        //- Reference to the points of sideA
        const pointField& pointsA_;

        //- Reference to the points of sideB
        const pointField& pointsB_;

        //- Flag to reverse B faces
        const bool reverseB_;


    // Static data members

        static scalar tol;


    // Private Member Functions

        //- Get triPoints from face
        inline triPoints getTriPoints
        (
            const pointField& points,
            const face& f,
            const bool reverse
        ) const;

        //- Set triPoints into tri list
        inline void setTriPoints
        (
            const point& a,
            const point& b,
            const point& c,
            label& count,
            FixedList<triPoints, 10>& tris
        ) const;

        //- Return point of intersection between plane and triangle edge
        inline point planeIntersection
        (
            const FixedList<scalar, 3>& d,
            const triPoints& t,
            const label negI,
            const label posI
        ) const;

        //- Return triangle area
        inline scalar triArea(const triPoints& t) const;


        //- Slice triangle with plane and generate new cut sub-triangles
        void triSliceWithPlane
        (
            const triPoints& tri,
            const plane& p,
            FixedList<triPoints, 10>& tris,
            label& nTris,
            const scalar len
        );

        //- Return area of intersection of triangles src and tgt
        scalar triangleIntersect
        (
            const triPoints& src,
            const triPoints& tgt,
            const vector& n
        );


public:

    // Constructors

        //- Construct from components
        faceAreaIntersect
        (
            const pointField& pointsA,
            const pointField& pointsB,
            const bool reverseB = false
        );


    // Public Member Functions

        //- Fraction of local length scale to use as intersection tolerance
        inline static scalar& tolerance();

        //- Triangulate a face using the given triangulation mode
        static void triangulate
        (
            const face& f,
            const pointField& points,
            const triangulationMode& triMode,
            faceList& faceTris
        );

        //- Return area of intersection of faceA with faceB
        scalar calc
        (
            const face& faceA,
            const face& faceB,
            const vector& n,
            const triangulationMode& triMode
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faceAreaIntersectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
