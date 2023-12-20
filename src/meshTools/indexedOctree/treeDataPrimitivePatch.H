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
    Foam::treeDataPrimitivePatch

Description
    Encapsulation of data needed to search on PrimitivePatches

SourceFiles
    treeDataPrimitivePatch.C

\*---------------------------------------------------------------------------*/

#ifndef treeDataPrimitivePatch_H
#define treeDataPrimitivePatch_H

#include "treeBoundBoxList.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class Type> class indexedOctree;


/*---------------------------------------------------------------------------*\
                 Class treeDataPrimitivePatchName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(treeDataPrimitivePatch);


/*---------------------------------------------------------------------------*\
                   Class treeDataPrimitivePatch Declaration
\*---------------------------------------------------------------------------*/

template<class PatchType>
class treeDataPrimitivePatch
:
    public treeDataPrimitivePatchName
{
    // Private Data

        //- Underlying geometry
        const PatchType& patch_;

        //- Whether to precalculate and store face bounding box
        const bool cacheBb_;

        //- Tolerance to use for intersection tests
        const scalar planarTol_;

        //- Face bounding boxes (valid only if cacheBb_)
        treeBoundBoxList bbs_;


    // Private Member Functions

        //- Calculate face bounding box
        static treeBoundBox calcBb(const pointField&, const face&);

        //- Initialise all member data
        void update();


public:


    class findNearestOp
    {
        const indexedOctree<treeDataPrimitivePatch>& tree_;

    public:

        findNearestOp(const indexedOctree<treeDataPrimitivePatch>& tree);

        void operator()
        (
            const labelUList& indices,
            const point& sample,

            scalar& nearestDistSqr,
            label& minIndex,
            point& nearestPoint
        ) const;

        //- Calculates nearest (to line) point in shape.
        //  Returns point and distance (squared)
        void operator()
        (
            const labelUList& indices,
            const linePointRef& ln,

            treeBoundBox& tightest,
            label& minIndex,
            point& linePoint,
            point& nearestPoint
        ) const;
    };


    class findIntersectOp
    {
        const indexedOctree<treeDataPrimitivePatch>& tree_;

    public:

        findIntersectOp(const indexedOctree<treeDataPrimitivePatch>& tree);

        //- Calculate intersection of any face with ray. Sets result
        //  accordingly. Used to find first intersection.
        bool operator()
        (
            const label index,
            const point& start,
            const point& end,
            point& intersectionPoint
        ) const;
    };


    class findAllIntersectOp
    {
        const indexedOctree<treeDataPrimitivePatch>& tree_;

        DynamicList<label>& shapeMask_;

    public:

        findAllIntersectOp
        (
            const indexedOctree<treeDataPrimitivePatch>& tree,
            DynamicList<label>& shapeMask
        );

        //- Calculate intersection of unique face with ray. Sets result
        //  accordingly. Used to find all faces.
        bool operator()
        (
            const label index,
            const point& start,
            const point& end,
            point& intersectionPoint
        ) const;
    };


    class findSelfIntersectOp
    {
        const indexedOctree<treeDataPrimitivePatch>& tree_;

        const label edgeIndex_;

    public:

        findSelfIntersectOp
        (
            const indexedOctree<treeDataPrimitivePatch>& tree,
            const label edgeID
        );

        //- Calculate intersection of face with edge of patch. Excludes
        //  faces that use edgeID. Used to find self intersection.
        bool operator()
        (
            const label index,
            const point& start,
            const point& end,
            point& intersectionPoint
        ) const;
    };


    // Constructors

        //- Construct from patch.
        treeDataPrimitivePatch
        (
            const bool cacheBb,
            const PatchType&,
            const scalar planarTol
        );


    // Member Functions

        // Access

            label size() const
            {
                return patch_.size();
            }

            //- Get representative point cloud for all shapes inside
            //  (one point per shape)
            pointField shapePoints() const;

            //- Return access to the underlying patch
            const PatchType& patch() const
            {
                return patch_;
            }


        // Search

            //- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
            //  Only makes sense for closed surfaces.
            volumeType getVolumeType
            (
                const indexedOctree<treeDataPrimitivePatch<PatchType>>&,
                const point&
            ) const;

            //- Does shape at index overlap bb
            bool overlaps
            (
                const label index,
                const treeBoundBox& sampleBb
            ) const;

            //- Does shape at index overlap sphere
            bool overlaps
            (
                const label index,
                const point& centre,
                const scalar radiusSqr
            ) const;

            //- Helper: find intersection of line with shapes
            static bool findIntersection
            (
                const indexedOctree<treeDataPrimitivePatch<PatchType>>& tree,
                const label index,
                const point& start,
                const point& end,
                point& intersectionPoint
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "treeDataPrimitivePatch.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
