/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::meshSearch

Description
    Various (local, not parallel) searches on polyMesh;
    uses (demand driven) octree to search.

SourceFiles
    meshSearch.C

\*---------------------------------------------------------------------------*/

#ifndef meshSearch_H
#define meshSearch_H

#include "pointIndexHit.H"
#include "pointField.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class treeDataCell;
class treeDataFace;
template<class Type> class indexedOctree;
class treeBoundBox;

/*---------------------------------------------------------------------------*\
                         Class meshSearch Declaration
\*---------------------------------------------------------------------------*/

class meshSearch
{
    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Whether to use cell decomposition for all geometric tests
        const polyMesh::cellDecomposition cellDecompMode_;

        //- Data bounding box
        mutable autoPtr<treeBoundBox> overallBbPtr_;

        //- Demand driven octrees
        mutable autoPtr<indexedOctree<treeDataFace>> boundaryTreePtr_;
        mutable autoPtr<indexedOctree<treeDataCell>> cellTreePtr_;


    // Private Member Functions

        //- Updates nearestI, nearestDistSqr from any closer ones.
        static bool findNearer
        (
            const point& sample,
            const pointField& points,
            label& nearestI,
            scalar& nearestDistSqr
        );

        //- Updates nearestI, nearestDistSqr from any selected closer ones.
        static bool findNearer
        (
            const point& sample,
            const pointField& points,
            const labelList& indices,
            label& nearestI,
            scalar& nearestDistSqr
        );


        // Cells

            //- Nearest cell centre using octree
            label findNearestCellTree(const point&) const;

            //- Nearest cell centre going through all cells
            label findNearestCellLinear(const point&) const;

            //- Walk from seed. Does not 'go around' boundary, just returns
            //  last cell before boundary.
            label findNearestCellWalk(const point&, const label) const;

            //- Cell containing location. Linear search.
            label findCellLinear(const point&) const;

            //- Walk from seed. Does not 'go around' boundary, just returns
            //  last cell before boundary.
            label findCellWalk(const point&, const label) const;


        // Faces

            label findNearestFaceTree(const point&) const;

            label findNearestFaceLinear(const point&) const;

            label findNearestFaceWalk(const point&, const label) const;



        // Boundary faces

            //- Walk from seed to find nearest boundary face. Gets stuck in
            //  local minimum.
            label findNearestBoundaryFaceWalk
            (
                const point& location,
                const label seedFacei
            ) const;


public:

    // Declare name of the class and its debug switch
    ClassName("meshSearch");


    // Static Data Members

        //- Tolerance on linear dimensions
        static scalar tol_;


    // Constructors

        //- Construct from components. Constructs bb slightly bigger than
        //  mesh points bb.
        meshSearch
        (
            const polyMesh& mesh,
            const polyMesh::cellDecomposition = polyMesh::CELL_TETS
        );

        //- Construct with a custom bounding box. Any mesh element outside
        //  bb will not be found. Up to user to make sure bb
        //  extends slightly beyond wanted elements.
        meshSearch
        (
            const polyMesh& mesh,
            const treeBoundBox& bb,
            const polyMesh::cellDecomposition = polyMesh::CELL_TETS
        );

        //- Disallow default bitwise copy construction
        meshSearch(const meshSearch&) = delete;


    //- Destructor
    ~meshSearch();


    // Member Functions

        // Access

            const polyMesh& mesh() const
            {
                return mesh_;
            }

            polyMesh::cellDecomposition decompMode() const
            {
                return cellDecompMode_;
            }

            //- Get (demand driven) reference to octree holding all
            //  boundary faces
            const indexedOctree<treeDataFace>& boundaryTree() const;

            //- Get (demand driven) reference to octree holding all cells
            const indexedOctree<treeDataCell>& cellTree() const;


        // Queries

            //- Find nearest cell in terms of cell centre.
            //  Options:
            //  - use octree
            //  - use linear search
            //  - if seed is provided walk. (uses findNearestCellWalk;
            //    does not handle holes in domain)
            label findNearestCell
            (
                const point& location,
                const label seedCelli = -1,
                const bool useTreeSearch = true
            ) const;

            label findNearestFace
            (
                const point& location,
                const label seedFacei = -1,
                const bool useTreeSearch = true
            ) const;

            //- Find cell containing location.
            //  If seed provided walks and falls back to linear/tree search.
            //  (so handles holes correctly)s
            //  Returns -1 if not in domain.
            label findCell
            (
                const point& location,
                const label seedCelli = -1,
                const bool useTreeSearch = true
            ) const;

            //- Find nearest boundary face
            //  If seed provided walks but then does not pass local minima
            //  in distance. Also does not jump from one connected region to
            //  the next.
            label findNearestBoundaryFace
            (
                const point& location,
                const label seedFacei = -1,
                const bool useTreeSearch = true
            ) const;

            //- Find first intersection of boundary in segment [pStart, pEnd]
            //  (so inclusive of endpoints). Always octree for now
            pointIndexHit intersection(const point& pStart, const point& pEnd)
            const;

            //- Find all intersections of boundary within segment pStart .. pEnd
            //  Always octree for now
            List<pointIndexHit> intersections
            (
                const point& pStart,
                const point& pEnd
            ) const;

            //- Determine inside/outside status
            bool isInside(const point&) const;


        //- Delete all storage
        void clearOut();

        //- Correct for mesh geom/topo changes
        void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshSearch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
