/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::searchableSurfaces::plate

Description
    Surface geometry of a finite plane, aligned with the coordinate
    axes, which can be used with snappyHexMesh.

Usage
    \table
        Property  | Description                         | Required
        origin    | Position of one corner of the plate | yes
        span      | The vector span of the plate        | yes
    \endtable

    Note: the span of the plate requires one zero component corresponding to the
    normal direction of the plate.

    Example specification in snappyHexMeshDict/geometry:
    \verbatim
        type       plate;
        origin     (10 10 10);
        span       (20  0 10);
    \endverbatim

SourceFiles
    searchablePlate.C

\*---------------------------------------------------------------------------*/

#ifndef plate_searchableSurface_H
#define plate_searchableSurface_H

#include "searchableSurface.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace searchableSurfaces
{

/*---------------------------------------------------------------------------*\
                             Class plate Declaration
\*---------------------------------------------------------------------------*/

class plate
:
    public searchableSurface
{
    // Private Member Data

        const point origin_;

        const vector span_;

        //- Coordinate direction which is normal
        const direction normalDir_;

        mutable wordList regions_;


    // Private Member Functions

        //- Calculate normal direction from span
        static direction calcNormal(const point&);

        //- Inherit findNearest from searchableSurface
        using searchableSurface::findNearest;

        pointIndexHit findNearest
        (
            const point& sample,
            const scalar nearestDistSqr
        ) const;

        pointIndexHit findLine
        (
            const point& start,
            const point& end
        ) const;

        //- Disallow default bitwise copy construction
        plate(const plate&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const plate&) = delete;


public:

    //- Runtime type information
    TypeName("plate");


    // Constructors

        //- Construct from components
        plate
        (
            const IOobject& io,
            const point& origin,
            const vector& span
        );

        //- Construct from dictionary (used by searchableSurface)
        plate
        (
            const IOobject& io,
            const dictionary& dict
        );


    //- Destructor
    virtual ~plate();


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
        virtual tmp<pointField> coordinates() const;

        //- Get bounding spheres (centre and radius squared), one per element.
        //  Any point on element is guaranteed to be inside.
        virtual void boundingSpheres
        (
            pointField& centres,
            scalarField& radiusSqr
        ) const;

        //- Get the points that define the surface.
        virtual tmp<pointField> points() const;

        //- Does any part of the surface overlap the supplied bound box?
        virtual bool overlaps(const boundBox& bb) const;


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End of namespace searchableSurfaces
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
