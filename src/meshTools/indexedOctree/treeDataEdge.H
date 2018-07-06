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
    Foam::treeDataEdge

Description
    Holds data for octree to work on an edges subset.

SourceFiles
    treeDataEdge.C

\*---------------------------------------------------------------------------*/

#ifndef treeDataEdge_H
#define treeDataEdge_H

#include "treeBoundBoxList.H"
#include "linePointRef.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class Type> class indexedOctree;

/*---------------------------------------------------------------------------*\
                           Class treeDataEdge Declaration
\*---------------------------------------------------------------------------*/

class treeDataEdge
{
    // Static data

        //- Tolerance on linear dimensions
        static scalar tol;


    // Private data

        //- Reference to edgeList
        const edgeList& edges_;

        //- Reference to points
        const pointField& points_;

        //- Labels of edges
        const labelList edgeLabels_;

        //- Whether to precalculate and store face bounding box
        const bool cacheBb_;

        //- Bbs for all above edges (valid only if cacheBb_)
        treeBoundBoxList bbs_;


    // Private Member Functions

        //- Calculate edge bounding box
        treeBoundBox calcBb(const label edgeI) const;

        //- Initialise all member data
        void update();

public:


    class findNearestOp
    {
        const indexedOctree<treeDataEdge>& tree_;

    public:

        findNearestOp(const indexedOctree<treeDataEdge>& tree);

        void operator()
        (
            const labelUList& indices,
            const point& sample,

            scalar& nearestDistSqr,
            label& minIndex,
            point& nearestPoint
        ) const;

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
    public:

        findIntersectOp(const indexedOctree<treeDataEdge>& tree);

        //- Calculate intersection of triangle with ray. Sets result
        //  accordingly
        bool operator()
        (
            const label index,
            const point& start,
            const point& end,
            point& intersectionPoint
        ) const;
    };


    // Declare name of the class and its debug switch
    ClassName("treeDataEdge");


    // Constructors

        //- Construct from selected edges. !Holds references to edges and points
        treeDataEdge
        (
            const bool cacheBb,
            const edgeList& edges,
            const pointField& points,
            const labelUList& edgeLabels
        );

        //- Construct from selected edges, transferring contents.
        // !Holds references to edges and points
        treeDataEdge
        (
            const bool cacheBb,
            const edgeList& edges,
            const pointField& points,
            const Xfer<labelList>& edgeLabels
        );


    // Member Functions

        // Access

            const edgeList& edges() const
            {
                return edges_;
            }

            const pointField& points() const
            {
                return points_;
            }

            const labelList& edgeLabels() const
            {
                return edgeLabels_;
            }

            label size() const
            {
                return edgeLabels_.size();
            }

            //- Get representative point cloud for all shapes inside
            //  (one point per shape)
            pointField shapePoints() const;


        // Search

            //- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
            //  Only makes sense for closed surfaces.
            volumeType getVolumeType
            (
                const indexedOctree<treeDataEdge>&,
                const point&
            ) const;

            //- Does (bb of) shape at index overlap bb
            bool overlaps
            (
                const label index,
                const treeBoundBox& sampleBb
            ) const;

            //- Does (bb of) shape at index overlap bb
            bool overlaps
            (
                const label index,
                const point& centre,
                const scalar radiusSqr
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
