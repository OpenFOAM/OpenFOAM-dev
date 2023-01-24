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
    Foam::searchableSurfaceCollection

Description
    Makes a collection of surface geometries by copying from an existing
    defined surface geometry. There are no boolean operations, e.g.
    overlapping surfaces are not intersected.

Usage
    searchableSurfaceCollection requires a surface geometry to be defined
    initially, e.g. \c buildingA in the example below. The collection then
    defines copies of buildingA which each can be scaled using the mandatory
    \c scale parameter and then rotated and translated by the mandatory
    \c transform.  In the example below, two geometries are included named
    \c buildingB and \c buildingC which are both formed by a translation
    of \c buildingA according to the \c origin parameter. No rotation is
    applied (by setting e1 and e2 to the global x and y axis directions,
    respectively).

    Example specification in snappyHexMeshDict/geometry:

    \verbatim
    buildingA
    {
        type       searchableBox;
        min        (100 100  0);
        max        (120 120  100);
    }

    moreBuildings
    {
        type searchableSurfaceCollection;

        mergeSubRegions true;

        buildingB
        {
            surface buildingA;
            scale   (1 1 1);
            transform
            {
                coordinateSystem
                {
                    type cartesian;
                    coordinateRotation
                    {
                        type    axesRotation;
                        e1      (1 0 0);
                        e2      (0 1 0);
                    }
                    origin  (40 40 0);
                }
            }
        }

        buildingC
        {
            $buildingB;
            transform
            {
                coordinateSystem { origin  (30 -40 0); }
            }
        }
    }
    \endverbatim

SourceFiles
    searchableSurfaceCollection.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurfaceCollection_H
#define searchableSurfaceCollection_H

#include "searchableSurface.H"
#include "treeBoundBox.H"
#include "coordinateSystem.H"
#include "UPtrList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class searchableSurfaceCollection Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaceCollection
:
    public searchableSurface
{
    // Private Member Data

        // Per instance data

            //- Instance name
            wordList instance_;

            //- Scaling vector
            vectorField scale_;

            //- transformation
            PtrList<coordinateSystem> transform_;

            UPtrList<searchableSurface> subGeom_;

            Switch mergeSubRegions_;

            //- Offsets for indices coming from different surfaces
            //  (sized with size() of each surface)
            labelList indexOffset_;

        //- Region names
        mutable wordList regions_;

        //- From individual regions to collection regions
        mutable labelList regionOffset_;


    // Private Member Functions

        //- Inherit findNearest from searchableSurface
        using searchableSurface::findNearest;

        //- Find point nearest to sample. Updates minDistSqr. Sets nearestInfo
        //  and surface index
        void findNearest
        (
            const pointField& samples,
            scalarField& minDistSqr,
            List<pointIndexHit>& nearestInfo,
            labelList& nearestSurf
        ) const;

        //- Sort hits into per-surface bins. Misses are rejected.
        //  Maintains map back to position
        void sortHits
        (
            const List<pointIndexHit>& info,
            List<List<pointIndexHit>>& surfInfo,
            labelListList& infoMap
        ) const;


public:

    //- Runtime type information
    TypeName("searchableSurfaceCollection");


    // Constructors

        //- Construct from dictionary (used by searchableSurface)
        searchableSurfaceCollection
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        searchableSurfaceCollection
        (
            const searchableSurfaceCollection&
        ) = delete;


    //- Destructor
    virtual ~searchableSurfaceCollection();


    // Member Functions

        //- Scaling vector per subsurface
        const vectorField& scale() const
        {
            return scale_;
        }

        //- Scaling vector per subsurface
        vectorField& scale()
        {
            return scale_;
        }

        //- Coordinate system per subsurface
        const PtrList<coordinateSystem>& transform() const
        {
            return transform_;
        }

        //- Coordinate system per subsurface
        PtrList<coordinateSystem>& transform()
        {
            return transform_;
        }

        virtual const wordList& regions() const;

        //- Whether supports volume type below
        virtual bool hasVolumeType() const
        {
            return false;
        }

        //- Range of local indices that can be returned.
        virtual label size() const;

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

        // Other

            //- Set bounds of surface. Bounds currently set as list of
            //  bounding boxes. The bounds are hints to the surface as for
            //  the range of queries it can expect. faceMap/pointMap can be
            //  set if the surface has done any redistribution.
            virtual void distribute
            (
                const List<treeBoundBox>&,
                const bool keepNonLocal,
                autoPtr<distributionMap>& faceMap,
                autoPtr<distributionMap>& pointMap
            );

            //- WIP. Store element-wise field.
            virtual void setField(const labelList& values);

            //- WIP. From a set of hits (points and
            //  indices) get the specified field. Misses do not get set. Return
            //  empty field if not supported.
            virtual void getField(const List<pointIndexHit>&, labelList&) const;

        // regIOobject implementation

            bool writeData(Ostream&) const
            {
                NotImplemented;
                return false;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const searchableSurfaceCollection&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
