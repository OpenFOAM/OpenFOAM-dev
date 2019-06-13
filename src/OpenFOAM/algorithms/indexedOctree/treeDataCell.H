/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::treeDataCell

Description
    Encapsulation of data needed to search in/for cells. Used to find the
    cell containing a point (e.g. cell-cell mapping).

SourceFiles
    treeDataCell.C

\*---------------------------------------------------------------------------*/

#ifndef treeDataCell_H
#define treeDataCell_H

#include "polyMesh.H"
#include "treeBoundBoxList.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class Type> class indexedOctree;

/*---------------------------------------------------------------------------*\
                        Class treeDataCell Declaration
\*---------------------------------------------------------------------------*/

class treeDataCell
{
    // Private Data

        const polyMesh& mesh_;

        //- Subset of cells to work on
        const labelList cellLabels_;

        //- Whether to precalculate and store cell bounding box
        const bool cacheBb_;

        //- How to decide if point is inside cell
        const polyMesh::cellDecomposition decompMode_;

        //- Cell bounding boxes (valid only if cacheBb_)
        treeBoundBoxList bbs_;


    // Private Member Functions

        //- Calculate cell bounding box
        treeBoundBox calcCellBb(const label celli) const;

        //- Initialise all member data
        void update();

public:


    class findNearestOp
    {
        const indexedOctree<treeDataCell>& tree_;

    public:

        findNearestOp(const indexedOctree<treeDataCell>& tree);

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
        const indexedOctree<treeDataCell>& tree_;

    public:

        findIntersectOp(const indexedOctree<treeDataCell>& tree);

        bool operator()
        (
            const label index,
            const point& start,
            const point& end,
            point& intersectionPoint
        ) const;
    };


    // Declare name of the class and its debug switch
    ClassName("treeDataCell");


    // Constructors

        //- Construct from mesh and subset of cells.
        treeDataCell
        (
            const bool cacheBb,
            const polyMesh&,
            const labelUList&,
            const polyMesh::cellDecomposition decompMode
        );

        //- Move construct from mesh and subset of cells, transferring contents
        treeDataCell
        (
            const bool cacheBb,
            const polyMesh&,
            labelList&&,
            const polyMesh::cellDecomposition decompMode
        );

        //- Construct from mesh. Uses all cells in mesh.
        treeDataCell
        (
            const bool cacheBb,
            const polyMesh&,
            const polyMesh::cellDecomposition decompMode
        );


    // Member Functions

        // Access

            inline const labelList& cellLabels() const
            {
                return cellLabels_;
            }

            inline const polyMesh& mesh() const
            {
                return mesh_;
            }

            inline polyMesh::cellDecomposition decompMode() const
            {
                return decompMode_;
            }

            inline label size() const
            {
                return cellLabels_.size();
            }

            //- Get representative point cloud for all shapes inside
            //  (one point per shape)
            pointField shapePoints() const;


        // Search

            //- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
            //  Only makes sense for closed surfaces.
            volumeType getVolumeType
            (
                const indexedOctree<treeDataCell>&,
                const point&
            ) const
            {
                NotImplemented;
                return volumeType::unknown;
            }

            //- Does (bb of) shape at index overlap bb
            bool overlaps
            (
                const label index,
                const treeBoundBox& sampleBb
            ) const;

            //- Does shape at index contain sample
            bool contains
            (
                const label index,
                const point& sample
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
