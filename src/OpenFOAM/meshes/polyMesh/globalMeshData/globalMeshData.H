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
    Foam::globalMeshData

Description
    Various mesh related information for a parallel run. Upon construction,
    constructs all info using parallel communication.

    Requires:
    - all processor patches to have correct ordering.
    - all processorPatches to have their transforms set.

    The shared point and edge addressing calculates addressing for points
    and edges on coupled patches.  In the 'old' way a distinction was made
    between points/edges that are only on two processors and those that are
    on multiple processors.  The problem is that those on multiple processors
    do not allow any transformations and require a global reduction on the
    master processor.

    The alternative is to have an exchange schedule
    (through a 'distributionMap')
    which sends all point/edge data (no distinction is made between
    those on two and those on more than two coupled patches) to the local
    'master'.  This master then does any calculation and sends
    the result back to the 'slave' points/edges.  This only needs to be done
    on points on coupled faces.  Any transformation is done using a
    predetermined set of transformations - since transformations have to be
    space filling only a certain number of transformation is supported.

    The exchange needs
    - a field of data
    - a distributionMap which does all parallel exchange and transformations
      This appends remote data to the end of the field.
    - a set of indices which indicate where to get untransformed data in the
      field
    - a set of indices which indicate where to get transformed data in the
      field

    Note:
      - compared to 17x nTotalFaces, nTotalPoints do not compensate for
        shared points since this would trigger full connectivity analysis
      - most calculation is demand driven and uses parallel communication
        so make sure to invoke on all processors at the same time
      - old sharedEdge calculation: currently an edge is considered shared if it
        uses two shared points and is used more than once.  This is not correct
        on processor patches but it only slightly overestimates the number of
        shared edges.  Doing full analysis of how many patches use the edge
        would be too complicated

See also
    distributionMap
    globalIndexAndTransform

SourceFiles
    globalMeshData.C
    globalMeshDataTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef globalMeshData_H
#define globalMeshData_H

#include "processorTopology.H"
#include "labelPair.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class polyMesh;
class distributionMap;
template<class T> class EdgeMap;
class globalIndex;
class globalIndexAndTransform;
class PackedBoolList;

/*---------------------------------------------------------------------------*\
                       Class globalMeshData Declaration
\*---------------------------------------------------------------------------*/

class globalMeshData
{

    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;


        // Data related to the complete mesh

            //- Total number of points in the complete mesh
            label nTotalPoints_;

            //- Total number of faces in the complete mesh
            label nTotalFaces_;

            //- Total number of cells in the complete mesh
            label nTotalCells_;


        // Processor patch addressing (be careful if not running in parallel!)

            //- Processor topology class
            processorTopology processorTopology_;

            //- List of processor patch labels
            //  (size of list = number of processor patches)
            labelList processorPatches_;

            //- List of indices into processorPatches_ for each patch.
            //  Index = -1 for non-processor patches.
            //  (size of list = number of patches)
            labelList processorPatchIndices_;

            //- processorPatchIndices_ of the neighbours processor patches
            labelList processorPatchNeighbours_;


        // Coupled point addressing
        // This is addressing from coupled point to coupled points/faces/cells.
        // This is a full schedule so includes points used by only two
        // coupled patches.

            //- Patch of coupled faces. Additional patch edge to mesh edges
            //  correspondence:
            //      points: meshPoints(), meshPointMap()
            //      edges : meshEdges(), meshEdgeMap()
            mutable autoPtr<indirectPrimitivePatch> coupledPatchPtr_;
            mutable autoPtr<labelList> coupledPatchMeshEdgesPtr_;
            mutable autoPtr<Map<label>> coupledPatchMeshEdgeMapPtr_;

            //- Global numbering for coupledPatch points
            mutable autoPtr<globalIndex> globalPointNumberingPtr_;

            //- Global numbering for transforms
            mutable autoPtr<globalIndexAndTransform> globalTransformsPtr_;

            // Coupled point to coupled points

            mutable autoPtr<labelListList> globalPointSlavesPtr_;
            mutable autoPtr<labelListList> globalPointTransformedSlavesPtr_;
            mutable autoPtr<distributionMap> globalPointSlavesMapPtr_;

            // Coupled edge to coupled edges

            mutable autoPtr<globalIndex> globalEdgeNumberingPtr_;
            mutable autoPtr<labelListList> globalEdgeSlavesPtr_;
            mutable autoPtr<labelListList> globalEdgeTransformedSlavesPtr_;
            mutable autoPtr<PackedBoolList> globalEdgeOrientationPtr_;
            mutable autoPtr<distributionMap> globalEdgeSlavesMapPtr_;


            // Coupled point to boundary faces

            mutable autoPtr<globalIndex> globalBoundaryFaceNumberingPtr_;
            mutable autoPtr<labelListList> globalPointBoundaryFacesPtr_;
            mutable autoPtr<labelListList>
                globalPointTransformedBoundaryFacesPtr_;
            mutable autoPtr<distributionMap> globalPointBoundaryFacesMapPtr_;

            // Coupled point to boundary cells

            mutable autoPtr<labelList> boundaryCellsPtr_;
            mutable autoPtr<globalIndex> globalBoundaryCellNumberingPtr_;
            mutable autoPtr<labelListList> globalPointBoundaryCellsPtr_;
            mutable autoPtr<labelListList>
                globalPointTransformedBoundaryCellsPtr_;
            mutable autoPtr<distributionMap> globalPointBoundaryCellsMapPtr_;


            // Other: coupled point to coupled COLLOCATED points
            mutable autoPtr<labelListList> globalCoPointSlavesPtr_;
            mutable autoPtr<distributionMap> globalCoPointSlavesMapPtr_;



        // Globally shared point addressing

            //- Total number of global points
            mutable label nGlobalPoints_;

            //- Indices of local points that are globally shared
            mutable autoPtr<labelList> sharedPointLabelsPtr_;

            //- Indices of globally shared points in the master list
            //  This list contains all the shared points in the mesh
            mutable autoPtr<labelList> sharedPointAddrPtr_;

            //- Shared point global labels.
            //  Global point index for every local shared point.
            //  Only valid if constructed with this information or if
            //  pointProcAddressing read.
            mutable autoPtr<labelList> sharedPointGlobalLabelsPtr_;


        // Globally shared edge addressing. Derived from shared points.
        // All demand driven since don't want to construct edges always.

            //- Total number of global edges
            mutable label nGlobalEdges_;

            //- Indices of local edges that are globally shared
            mutable autoPtr<labelList> sharedEdgeLabelsPtr_;

            //- Indices of globally shared edge in the master list
            //  This list contains all the shared edges in the mesh
            mutable autoPtr<labelList> sharedEdgeAddrPtr_;


    // Private Member Functions

        //- Set up processor patch addressing
        void initProcAddr();

        //- Helper function for shared edge addressing
        static void countSharedEdges
        (
            const EdgeMap<labelList>&,
            EdgeMap<label>&,
            label&
        );

        //- Calculate shared point addressing
        void calcSharedPoints() const;

        //- Calculate shared edge addressing
        void calcSharedEdges() const;

        //- Calculate global point addressing.
        void calcGlobalPointSlaves() const;


        // Global edge addressing

            //- Calculate connected points
            void calcPointConnectivity(List<labelPairList>&) const;

            //- Calculate pointEdges and pointPoints addressing
            void calcGlobalPointEdges
            (
                labelListList& globalPointEdges,
                List<labelPairList>& globalPointPoints
            ) const;

            //- Look up remote and local point and find using info the
            //  transforms to go from remotePoint to localPoint
            label findTransform
            (
                const labelPairList& info,
                const labelPair& remotePoint,
                const label localPoint
            ) const;

            //- Calculate global edge addressing.
            void calcGlobalEdgeSlaves() const;

            //- Calculate orientation w.r.t. edge master.
            void calcGlobalEdgeOrientation() const;


        // Global boundary face/cell addressing

            //- Calculate coupled point to uncoupled boundary faces. Local only.
            void calcPointBoundaryFaces(labelListList&) const;

            //- Calculate global point to global boundary face addressing.
            void calcGlobalPointBoundaryFaces() const;

            //- Calculate global point to global boundary cell addressing.
            void calcGlobalPointBoundaryCells() const;


        // Other

            // Point to collocated points. Note that not all points on
            // coupled patches now have a master! (since points on either
            // side of a cyclic are not connected). So check whether the map
            // reaches all points and decide who is master, slave and who is
            // its own master. Maybe store as well?

            void calcGlobalCoPointSlaves() const;


public:

    // Public class

        // To combineReduce a List. Just appends all lists.
        template<class T>
        class ListPlusEqOp
        {

        public:

            void operator()(T& x, const T& y) const
            {
                label n = x.size();

                x.setSize(x.size() + y.size());

                forAll(y, i)
                {
                    x[n++] = y[i];
                }
            }
        };


    //- Runtime type information
    ClassName("globalMeshData");


    // Static Data Members

        //- Geometric tolerance (fraction of bounding box)
        static const Foam::scalar matchTol_;


    // Constructors

        //- Construct from mesh, derive rest (does parallel communication!)
        globalMeshData(const polyMesh& mesh);

        //- Disallow default bitwise copy construction
        globalMeshData(const globalMeshData&) = delete;


    //- Destructor
    ~globalMeshData();

        //- Remove all demand driven data
        void clearOut();


    // Member Functions

        // Access

            //- Return the mesh reference
            const polyMesh& mesh() const
            {
                return mesh_;
            }

            //- Does the mesh contain processor patches? (also valid when
            //  not running parallel)
            bool parallel() const
            {
                return processorPatches_.size() > 0;
            }

            //- Return total number of points in decomposed mesh. Not
            //  compensated for duplicate points!
            label nTotalPoints() const
            {
                return nTotalPoints_;
            }

            //- Return total number of faces in decomposed mesh. Not
            //  compensated for duplicate faces!
            label nTotalFaces() const
            {
                return nTotalFaces_;
            }

            //- Return total number of cells in decomposed mesh.
            label nTotalCells() const
            {
                return nTotalCells_;
            }


        // Processor patch addressing (be careful when not running in parallel)

            //- Return the processor-processor connection table
            const labelListList& procNbrProcs() const
            {
                return processorTopology_.procNbrProcs();
            }

            //- Order in which the patches should be initialised/evaluated
            //  corresponding to the schedule
            const lduSchedule& patchSchedule() const
            {
                return processorTopology_.patchSchedule();
            }

            //- Return list of processor patch labels
            //  (size of list = number of processor patches)
            const labelList& processorPatches() const
            {
                return processorPatches_;
            }

            //- Return list of indices into processorPatches_ for each patch.
            //  Index = -1 for non-processor parches.
            //  (size of list = number of patches)
            const labelList& processorPatchIndices() const
            {
                return processorPatchIndices_;
            }

            //- Return processorPatchIndices of the neighbours
            //  processor patches. -1 if not running parallel.
            const labelList& processorPatchNeighbours() const
            {
                return processorPatchNeighbours_;
            }


        // Globally shared point addressing

            //- Return number of globally shared points
            label nGlobalPoints() const;

            //- Return indices of local points that are globally shared
            const labelList& sharedPointLabels() const;

            //- Return addressing into the complete globally shared points
            //  list
            //  Note: It is assumed that a (never constructed) complete
            //  list of globally shared points exists.  The set of shared
            //  points on the current processor is a subset of all shared
            //  points. Shared point addressing gives the index in the
            //  list of all globally shared points for each of the locally
            //  shared points.
            const labelList& sharedPointAddr() const;

            //- Return shared point global labels. Tries to read
            //  'pointProcAddressing' and returns list or -1 if none
            //  available.
            const labelList& sharedPointGlobalLabels() const;

            //- Collect coordinates of shared points on all processors.
            //  Does parallel communication. Not valid for points on
            //  cyclics since these combine into a single shared point but
            //  span multiple locations.
            pointField sharedPoints() const;



        // Globally shared edge addressing

            //- Return number of globally shared edges. Demand-driven
            //  calculation so call needs to be synchronous among processors!
            label nGlobalEdges() const;

            //- Return indices of local edges that are globally shared.
            //  Demand-driven
            //  calculation so call needs to be synchronous among processors!
            const labelList& sharedEdgeLabels() const;

            //- Return addressing into the complete globally shared edge
            //  list. The set of shared
            //  edges on the current processor is a subset of all shared
            //  edges. Shared edge addressing gives the index in the
            //  list of all globally shared edges for each of the locally
            //  shared edges.
            //  Demand-driven
            //  calculation so call needs to be synchronous among processors!
            const labelList& sharedEdgeAddr() const;



        // Global master - slave point communication

            //- Return patch of all coupled faces
            const indirectPrimitivePatch& coupledPatch() const;

            //- Return map from coupledPatch edges to mesh edges
            const labelList& coupledPatchMeshEdges() const;

            //- Return map from mesh edges to coupledPatch edges
            const Map<label>& coupledPatchMeshEdgeMap() const;

            //- Global transforms numbering
            const globalIndexAndTransform& globalTransforms() const;

            //- Helper: synchronise data with transforms
            template<class Type, class CombineOp, class TransformOp>
            static void syncData
            (
                List<Type>& pointData,
                const labelListList& slaves,
                const labelListList& transformedSlaves,
                const distributionMap& slavesMap,
                const globalIndexAndTransform&,
                const CombineOp& cop,
                const TransformOp& top
            );

            //- Helper: synchronise data without transforms
            template<class Type, class CombineOp>
            static void syncData
            (
                List<Type>& pointData,
                const labelListList& slaves,
                const labelListList& transformedSlaves,
                const distributionMap& slavesMap,
                const CombineOp& cop
            );


            // Coupled point to coupled points. Coupled points are
            // points on any coupled patch.

                //- Numbering of coupled points is according to coupledPatch.
                const globalIndex& globalPointNumbering() const;
                const labelListList& globalPointSlaves() const;
                const labelListList& globalPointTransformedSlaves() const;
                const distributionMap& globalPointSlavesMap() const;
                //- Helper to synchronise coupled patch point data
                template<class Type, class CombineOp, class TransformOp>
                void syncPointData
                (
                    List<Type>& pointData,
                    const CombineOp& cop,
                    const TransformOp& top
                ) const;

            // Coupled edge to coupled edges.

                const globalIndex& globalEdgeNumbering() const;
                const labelListList& globalEdgeSlaves() const;
                const labelListList& globalEdgeTransformedSlaves() const;
                const distributionMap& globalEdgeSlavesMap() const;
                //- Is my edge same orientation as master edge
                const PackedBoolList& globalEdgeOrientation() const;

            // Collocated point to collocated point

                const labelListList& globalCoPointSlaves() const;
                const distributionMap& globalCoPointSlavesMap() const;

            // Coupled point to boundary faces. These are uncoupled boundary
            // faces only but include empty patches.

                //- Numbering of boundary faces is face-mesh.nInternalFaces()
                const globalIndex& globalBoundaryFaceNumbering() const;
                const labelListList& globalPointBoundaryFaces() const;
                const labelListList& globalPointTransformedBoundaryFaces()
                const;
                const distributionMap& globalPointBoundaryFacesMap() const;

            // Coupled point to boundary cell

                //- From boundary cell to mesh cell
                const labelList& boundaryCells() const;

                //- Numbering of boundary cells is according to boundaryCells()
                const globalIndex& globalBoundaryCellNumbering() const;
                const labelListList& globalPointBoundaryCells() const;
                const labelListList& globalPointTransformedBoundaryCells()
                const;
                const distributionMap& globalPointBoundaryCellsMap() const;


            // Other

                //- Helper for merging (collocated!) mesh point data.
                //  Determines:
                //  - my unique indices
                //  - global numbering over all unique indices
                //  - the global number for all local points (so this will
                //    be local for my unique points)
                autoPtr<globalIndex> mergePoints
                (
                    labelList& pointToGlobal,
                    labelList& uniquePoints
                ) const;

                //- Helper for merging (collocated!) patch point data.
                //  Takes maps from:
                //  local points to/from mesh. Determines
                //  - my unique points. These are mesh point indices, not patch
                //    point indices.
                //  - global numbering over all unique indices.
                //  - the global number for all local points.
                autoPtr<globalIndex> mergePoints
                (
                    const labelList& meshPoints,
                    const Map<label>& meshPointMap,
                    labelList& pointToGlobal,
                    labelList& uniqueMeshPoints
                ) const;


        // Edit

            //- Update for moving points.
            void movePoints(const pointField& newPoints);

            //- Change global mesh data given a topological change. Does a
            //  full parallel analysis to determine shared points and
            //  boundaries.
            void topoChange();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const globalMeshData&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "globalMeshDataTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
