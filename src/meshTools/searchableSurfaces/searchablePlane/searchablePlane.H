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
    Foam::searchablePlane

Description
    Surface geometry of an infinite plane, which can be used with
    snappyHexMesh.

Usage
    \table
        Property  | Description                             | Required
        planeType | Plane type, e.g. pointAndNormal         | yes
        point     | Point location for pointAndNormal plane | yes
        normal    | Normal vector for pointAndNormal plane  | yes
    \endtable

    Note: there are also other options for planeType, including
    'planeEquation' and 'embeddedPoints'.

    Example specification in snappyHexMeshDict/geometry:
    \verbatim
        type       searchablePlane;
        planeType  pointAndNormal;
        point      (10 10 10);
        normal     (0  1  0);
    \endverbatim

SourceFiles
    searchablePlane.C

\*---------------------------------------------------------------------------*/

#ifndef searchablePlane_H
#define searchablePlane_H

#include "searchableSurface.H"
#include "plane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                       Class searchablePlane Declaration
\*---------------------------------------------------------------------------*/

class searchablePlane
:
    public searchableSurface,
    public plane
{
    // Private Member Data

        mutable wordList regions_;


    // Private Member Functions

        pointIndexHit findLine
        (
            const point& start,
            const point& end
        ) const;

        //- Return the boundBox of the plane
        boundBox calcBounds() const;


public:

    //- Runtime type information
    TypeName("searchablePlane");


    // Constructors

        //- Construct from components
        searchablePlane
        (
            const IOobject& io,
            const point& basePoint,
            const vector& normal
        );

        //- Construct from dictionary (used by searchableSurface)
        searchablePlane
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        searchablePlane(const searchablePlane&) = delete;


    //- Destructor
    virtual ~searchablePlane();


    // Member Functions

        virtual const wordList& regions() const;

        //- Whether supports volume type below
        virtual bool hasVolumeType() const
        {
            return false;
        }

        //- Range of local indices that can be returned.
        virtual label size() const
        {
            return 1;
        }

        //- Get representative set of element coordinates
        //  Usually the element centres (should be of length size()).
        virtual tmp<pointField> coordinates() const
        {
            tmp<pointField> tCtrs(new pointField(1, refPoint()));
            return tCtrs;
        }

        //- Get bounding spheres (centre and radius squared), one per element.
        //  Any point on element is guaranteed to be inside.
        //  Note: radius limited to sqr(great)
        virtual void boundingSpheres
        (
            pointField& centres,
            scalarField& radiusSqr
        ) const;

        //- Get the points that define the surface.
        virtual tmp<pointField> points() const
        {
            return coordinates();
        }

        //- Does any part of the surface overlap the supplied bound box?
        virtual bool overlaps(const boundBox& bb) const
        {
            NotImplemented;

            return false;
        }


        // Multiple point queries.

            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>&
            ) const;

            virtual void findLine
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            virtual void findLineAny
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            //- Get all intersections in order from start to end.
            virtual void findLineAll
            (
                const pointField& start,
                const pointField& end,
                List<List<pointIndexHit>>&
            ) const;

            //- From a set of points and indices get the region
            virtual void getRegion
            (
                const List<pointIndexHit>&,
                labelList& region
            ) const;

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Determine type (inside/outside/mixed) for point. unknown if
            //  cannot be determined (e.g. non-manifold surface)
            virtual void getVolumeType
            (
                const pointField&,
                List<volumeType>&
            ) const;


        // regIOobject implementation

            bool writeData(Ostream&) const
            {
                NotImplemented;
                return false;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const searchablePlane&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
