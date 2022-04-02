/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
class polyPatch;
class polyMesh;
class transformer;

// Forward declaration of friend functions and operators
template<class Derived>
class WallFaceBase;
template<class Derived>
Ostream& operator<<(Ostream&, const WallFaceBase<Derived>&);
template<class Derived>
Istream& operator>>(Istream&, WallFaceBase<Derived>&);

/*---------------------------------------------------------------------------*\
                          Class WallFaceBase Declaration
\*---------------------------------------------------------------------------*/

template<class Derived>
class WallFaceBase
{
    // Private Data

        //- Face points
        pointField points_;

        //- Normal distance (squared) from cell center to face
        scalar distSqr_;


protected:

    // Protected Member Functions

        //- ...
        template<class TrackingData>
        inline bool update
        (
            const point& pt,
            const WallFaceBase<Derived>& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    // Constructors

        //- Construct null
        inline WallFaceBase();

        //- Construct from face, distance
        inline WallFaceBase
        (
            const face& f,
            const pointField& points,
            const scalar distSqr
        );

        //- Construct from face, distance
        inline WallFaceBase
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

            //- Check whether the WallFaceBase has been changed at all or still
            //  contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data. Used for checking
            //  consistency across cyclics.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const polyMesh&,
                const WallFaceBase<Derived>&,
                const scalar,
                TrackingData& td
            ) const;

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const polyPatch& patch,
                const label patchFacei,
                const transformer& transform,
                TrackingData& td
            );

            //- Influence of neighbouring face
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh&,
                const label thisCelli,
                const label neighbourFacei,
                const WallFaceBase<Derived>& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const label neighbourCelli,
                const WallFaceBase<Derived>& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const WallFaceBase<Derived>& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Test equality
            template<class TrackingData>
            inline bool equal
            (
                const WallFaceBase<Derived>&,
                TrackingData& td
            ) const;


    // Member Operators

        inline bool operator==(const WallFaceBase<Derived>&) const;
        inline bool operator!=(const WallFaceBase<Derived>&) const;


    // IOstream Operators

        friend Ostream& operator<< <Derived>
        (
            Ostream&,
            const WallFaceBase<Derived>&
        );

        friend Istream& operator>> <Derived>
        (
            Istream&,
            WallFaceBase<Derived>&
        );
};

/*---------------------------------------------------------------------------*\
                          Class wallFace Declaration
\*---------------------------------------------------------------------------*/

class wallFace
:
    public WallFaceBase<wallFace>
{
    public:

        using WallFaceBase<wallFace>::WallFaceBase;

        template<class Derived> using type = WallFaceBase<Derived>;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallFaceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
