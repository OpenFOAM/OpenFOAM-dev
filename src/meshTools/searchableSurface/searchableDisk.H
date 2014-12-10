/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::searchableDisk

Description
    Searching on circular disk given as origin, normal (gets normalised)
    and radius

SourceFiles
    searchableDisk.C

\*---------------------------------------------------------------------------*/

#ifndef searchableDisk_H
#define searchableDisk_H

#include "treeBoundBox.H"
#include "searchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class searchableDisk Declaration
\*---------------------------------------------------------------------------*/

class searchableDisk
:
    public searchableSurface
{
private:

    // Private Member Data

        //- origin
        const point origin_;

        //- normal
        vector normal_;

        //- radius
        const scalar radius_;

        //- Names of regions
        mutable wordList regions_;


    // Private Member Functions

        //- Find nearest point on disk
        pointIndexHit findNearest
        (
            const point& sample,
            const scalar nearestDistSqr
        ) const;

        //- Find intersection with disk
        void findLine
        (
            const point& start,
            const point& end,
            pointIndexHit&
        ) const;

        //- Disallow default bitwise copy construct
        searchableDisk(const searchableDisk&);

        //- Disallow default bitwise assignment
        void operator=(const searchableDisk&);

public:

    //- Runtime type information
    TypeName("searchableDisk");


    // Constructors

        //- Construct from components
        searchableDisk
        (
            const IOobject& io,
            const point& origin,
            const point& normal,
            const scalar radius
        );

        //- Construct from dictionary (used by searchableSurface)
        searchableDisk
        (
            const IOobject& io,
            const dictionary& dict
        );

    //- Destructor
    virtual ~searchableDisk();


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
            tmp<pointField> tCtrs(new pointField(1, origin_));
            return tCtrs;
        }

        //- Get bounding spheres (centre and radius squared), one per element.
        //  Any point on element is guaranteed to be inside.
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
            notImplemented
            (
                "searchableDisk::overlaps(const boundBox&) const"
            );

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
                List<List<pointIndexHit> >&
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
                notImplemented("searchableDisk::writeData(Ostream&) const");
                return false;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
