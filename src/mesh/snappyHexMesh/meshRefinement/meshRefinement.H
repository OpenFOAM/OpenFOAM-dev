/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::meshRefinement

Description
    Helper class which maintains intersections of (changing) mesh with
    (static) surfaces.

    Maintains
    - per face any intersections of the cc-cc segment with any of the surfaces

SourceFiles
    meshRefinement.C
    meshRefinementBaffles.C
    meshRefinementMerge.C
    meshRefinementProblemCells.C
    meshRefinementRefine.C

\*---------------------------------------------------------------------------*/

#ifndef meshRefinement_H
#define meshRefinement_H

#include "hexRef8.H"
#include "mapPolyMesh.H"
#include "autoPtr.H"
#include "labelPair.H"
#include "indirectPrimitivePatch.H"
#include "pointFieldsFwd.H"
#include "Tuple2.H"
#include "pointIndexHit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class fvMesh;
class mapDistributePolyMesh;
class decompositionMethod;
class refinementSurfaces;
class refinementFeatures;
class refinementRegions;
class removeCells;
class fvMeshDistribute;
class searchableSurface;
class regionSplit;
class globalIndex;
class removePoints;
class localPointRegion;

class snapParameters;

/*---------------------------------------------------------------------------*\
                       Class meshRefinement Declaration
\*---------------------------------------------------------------------------*/

class meshRefinement
{

public:

    // Public data types

        //- Enumeration for what to debug
        enum IOdebugType
        {
            IOMESH,
            IOOBJINTERSECTIONS,
            IOFEATURESEEDS,
            IOATTRACTION,
            IOLAYERINFO
        };

        static const NamedEnum<IOdebugType, 5> IOdebugTypeNames;
        enum debugType
        {
            MESH = 1<<IOMESH,
            OBJINTERSECTIONS = 1<<IOOBJINTERSECTIONS,
            FEATURESEEDS = 1<<IOFEATURESEEDS,
            ATTRACTION = 1<< IOATTRACTION,
            LAYERINFO = 1<<IOLAYERINFO
        };

        //- Enumeration for what to output
        enum IOoutputType
        {
            IOOUTPUTLAYERINFO
        };

        static const NamedEnum<IOoutputType, 1> IOoutputTypeNames;
        enum outputType
        {
            OUTPUTLAYERINFO = 1<<IOOUTPUTLAYERINFO
        };

        //- Enumeration for what to write
        enum IOwriteType
        {
            IOWRITEMESH,
            IONOWRITEREFINEMENT,
            IOWRITELEVELS,
            IOWRITELAYERSETS,
            IOWRITELAYERFIELDS
        };

        static const NamedEnum<IOwriteType, 5> IOwriteTypeNames;
        enum writeType
        {
            WRITEMESH = 1<<IOWRITEMESH,
            NOWRITEREFINEMENT = 1<<IONOWRITEREFINEMENT,
            WRITELEVELS = 1<<IOWRITELEVELS,
            WRITELAYERSETS = 1<<IOWRITELAYERSETS,
            WRITELAYERFIELDS = 1<<IOWRITELAYERFIELDS
        };

        //- Enumeration for how the userdata is to be mapped upon refinement.
        enum mapType
        {
            MASTERONLY = 1, //!< maintain master only
            KEEPALL = 2,    //!< have slaves (upon refinement) from master
            REMOVE = 4      //!< set value to -1 any face that was refined
        };


private:

    // Static Data Members

        //- Control of writing level
        static writeType writeLevel_;

        //- Control of output/log level
        static outputType outputLevel_;


    // Private Data

        //- Reference to mesh
        fvMesh& mesh_;

        //- Tolerance used for sorting coordinates (used in 'less' routine)
        const scalar mergeDistance_;

        //- Overwrite the mesh?
        const bool overwrite_;

        //- Instance of mesh upon construction. Used when in overwrite_ mode.
        const word oldInstance_;

        //- All surface-intersection interaction
        const refinementSurfaces& surfaces_;

        //- All feature-edge interaction
        const refinementFeatures& features_;

        //- All shell-refinement interaction
        const refinementRegions& shells_;

        //- Refinement engine
        hexRef8 meshCutter_;

        //- Per cc-cc vector the index of the surface hit
        labelIOList surfaceIndex_;

        //- User supplied face based data.
        List<Tuple2<mapType, labelList>> userFaceData_;

        //- Meshed patches - are treated differently. Stored as wordList since
        //  order changes.
        wordList meshedPatches_;


    // Private Member Functions

        //- Find out which faces have changed given cells (old mesh labels)
        //  that were marked for refinement.
        static labelList getChangedFaces
        (
            const mapPolyMesh&,
            const labelList& oldCellsToRefine
        );

        //- Calculate coupled boundary end vector and refinement level
        void calcNeighbourData
        (
            labelList& neiLevel,
            pointField& neiCc
        ) const;

        //- Find any intersection of surface. Store in surfaceIndex_.
        void updateIntersections(const labelList& changedFaces);

        //- Remove cells. Put exposedFaces into exposedPatchIDs.
        autoPtr<mapPolyMesh> doRemoveCells
        (
            const labelList& cellsToRemove,
            const labelList& exposedFaces,
            const labelList& exposedPatchIDs,
            removeCells& cellRemover
        );

        // Get cells which are inside any closed surface. Note that
        // all closed surfaces
        // will have already been oriented to have keepPoint outside.
        labelList getInsideCells(const word&) const;

        // Do all to remove inside cells
        autoPtr<mapPolyMesh> removeInsideCells
        (
            const string& msg,
            const label exposedPatchi
        );


        // Refinement candidate selection

            //- Mark cell for refinement (if not already marked). Return false
            // if refinelimit hit. Keeps running count (in nRefine) of cells
            // marked for refinement
            static bool markForRefine
            (
                const label markValue,
                const label nAllowRefine,
                label& cellValue,
                label& nRefine
            );

            //- Mark every cell with level of feature passing through it
            //  (or -1 if not passed through). Uses tracking.
            void markFeatureCellLevel
            (
                const pointField& keepPoints,
                labelList& maxFeatureLevel
            ) const;

            //- Calculate list of cells to refine based on intersection of
            //  features.
            label markFeatureRefinement
            (
                const pointField& keepPoints,
                const label nAllowRefine,

                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark cells for distance-to-feature based refinement.
            label markInternalDistanceToFeatureRefinement
            (
                const label nAllowRefine,
                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark cells for refinement-shells based refinement.
            label markInternalRefinement
            (
                const label nAllowRefine,
                labelList& refineCell,
                label& nRefine
            ) const;

            //- Collect faces that are intersected and whose neighbours aren't
            //  yet marked  for refinement.
            labelList getRefineCandidateFaces
            (
                const labelList& refineCell
            ) const;

            //- Mark cells for surface intersection based refinement.
            label markSurfaceRefinement
            (
                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,
                labelList& refineCell,
                label& nRefine
            ) const;

            //- Helper: count number of normals1 that are in normals2
            label countMatches
            (
                const List<point>& normals1,
                const List<point>& normals2,
                const scalar tol = 1e-6
            ) const;

            //- Mark cells for surface curvature based refinement. Marks if
            //  local curvature > curvature or if on different regions
            //  (markDifferingRegions)
            label markSurfaceCurvatureRefinement
            (
                const scalar curvature,
                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,
                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark cell if local a gap topology or
            bool checkProximity
            (
                const scalar planarCos,
                const label nAllowRefine,

                const label surfaceLevel,
                const vector& surfaceLocation,
                const vector& surfaceNormal,

                const label celli,

                label& cellMaxLevel,
                vector& cellMaxLocation,
                vector& cellMaxNormal,

                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark cells for surface proximity based refinement.
            label markProximityRefinement
            (
                const scalar curvature,
                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,

                labelList& refineCell,
                label& nRefine
            ) const;

        // Baffle handling

            //- Get faces to repatch. Returns map from face to patch.
            Map<labelPair> getZoneBafflePatches
            (
                const bool allowBoundary,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch
            ) const;

            //- Determine patches for baffles on all intersected unnamed faces
            void getBafflePatches
            (
                const labelList& globalToMasterPatch,
                const labelList& neiLevel,
                const pointField& neiCc,
                labelList& ownPatch,
                labelList& nbrPatch
            ) const;

            //- Repatches external face or creates baffle for internal face
            //  with user specified patches (might be different for both sides).
            //  Returns label of added face.
            label createBaffle
            (
                const label facei,
                const label ownPatch,
                const label nbrPatch,
                polyTopoChange& meshMod
            ) const;

        // Problem cell handling

            //- Helper function to mark face as being on 'boundary'. Used by
            //  markFacesOnProblemCells
            void markBoundaryFace
            (
                const label facei,
                boolList& isBoundaryFace,
                boolList& isBoundaryEdge,
                boolList& isBoundaryPoint
            ) const;

            void findNearest
            (
                const labelList& meshFaces,
                List<pointIndexHit>& nearestInfo,
                labelList& nearestSurface,
                labelList& nearestRegion,
                vectorField& nearestNormal
            ) const;

            Map<label> findEdgeConnectedProblemCells
            (
                const scalarField& perpendicularAngle,
                const labelList&
            ) const;

            bool isCollapsedFace
            (
                const pointField&,
                const pointField& neiCc,
                const scalar minFaceArea,
                const scalar maxNonOrtho,
                const label facei
            ) const;

            bool isCollapsedCell
            (
                const pointField&,
                const scalar volFraction,
                const label celli
            ) const;

            //- Returns list with for every internal face -1 or the patch
            //  they should be baffled into. If removeEdgeConnectedCells is set
            //  removes cells based on perpendicularAngle.
            labelList markFacesOnProblemCells
            (
                const dictionary& motionDict,
                const bool removeEdgeConnectedCells,
                const scalarField& perpendicularAngle,
                const labelList& globalToMasterPatch
            ) const;

            //- Returns list with for every face the label of the nearest
            //  patch. Any unreached face (disconnected mesh?) becomes
            //  adaptPatchIDs[0]
            labelList nearestPatch(const labelList& adaptPatchIDs) const;

            //- Returns list with for every internal face -1 or the patch
            //  they should be baffled into.
            labelList markFacesOnProblemCellsGeometric
            (
                const snapParameters& snapParams,
                const dictionary& motionDict
            ) const;


        // Baffle merging

            //- Extract those baffles (duplicate) faces that are on the edge
            //  of a baffle region. These are candidates for merging.
            List<labelPair> freeStandingBaffles
            (
                const List<labelPair>&,
                const scalar freeStandingAngle
            ) const;


        // Zone handling

            //- Finds zone per cell for cells inside closed named surfaces.
            //  (uses geometric test for insideness)
            //  Adapts namedSurfaceIndex so all faces on boundary of cellZone
            //  have corresponding faceZone.
            void findCellZoneGeometric
            (
                const pointField& neiCc,
                const labelList& closedNamedSurfaces,
                labelList& namedSurfaceIndex,
                const labelList& surfaceToCellZone,
                labelList& cellToZone
            ) const;

            //- Finds zone per cell for cells inside named surfaces that have
            //  an inside point specified.
            void findCellZoneInsideWalk
            (
                const labelList& locationSurfaces,
                const labelList& namedSurfaceIndex,
                const labelList& surfaceToCellZone,
                labelList& cellToZone
            ) const;

            //- Determines cell zone from cell region information.
            bool calcRegionToZone
            (
                const label surfZoneI,
                const label ownRegion,
                const label neiRegion,

                labelList& regionToCellZone
            ) const;

            //- Finds zone per cell. Uses topological walk with all faces
            //  marked in namedSurfaceIndex regarded as blocked.
            void findCellZoneTopo
            (
                const point& keepPoint,
                const labelList& namedSurfaceIndex,
                const labelList& surfaceToCellZone,
                labelList& cellToZone
            ) const;

            //- Make namedSurfaceIndex consistent with cellToZone
            //  - clear out any blocked faces in between same cell zone.
            void makeConsistentFaceIndex
            (
                const labelList& cellToZone,
                labelList& namedSurfaceIndex
            ) const;

            //- Remove any loose standing cells
            void handleSnapProblems
            (
                const snapParameters& snapParams,
                const bool useTopologicalSnapDetection,
                const bool removeEdgeConnectedCells,
                const scalarField& perpendicularAngle,
                const dictionary& motionDict,
                Time& runTime,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch
            );


            // Some patch utilities

            //- Get all faces in faceToZone that have no cellZone on
            //  either side.
            labelList freeStandingBaffleFaces
            (
                const labelList& faceToZone,
                const labelList& cellToZone,
                const labelList& neiCellZone
            ) const;

            //- Determine per patch edge the number of master faces. Used
            //  to detect non-manifold situations.
            void calcPatchNumMasterFaces
            (
                const PackedBoolList& isMasterFace,
                const indirectPrimitivePatch& patch,
                labelList& nMasterFaces
            ) const;

            //- Determine per patch face the (singly-) connected zone it
            //  is in. Return overall number of zones.
            label markPatchZones
            (
                const indirectPrimitivePatch& patch,
                const labelList& nMasterFaces,
                labelList& faceToZone
            ) const;

            //- Make faces consistent.
            void consistentOrientation
            (
                const PackedBoolList& isMasterFace,
                const indirectPrimitivePatch& patch,
                const labelList& nMasterFaces,
                const labelList& faceToZone,
                const Map<label>& zoneToOrientation,
                boolList& meshFlipMap
            ) const;


public:

    //- Runtime type information
    ClassName("meshRefinement");


    // Constructors

        //- Construct from components
        meshRefinement
        (
            fvMesh& mesh,
            const scalar mergeDistance,
            const bool overwrite,
            const refinementSurfaces&,
            const refinementFeatures&,
            const refinementRegions&
        );

        //- Disallow default bitwise copy construction
        meshRefinement(const meshRefinement&) = delete;


    // Member Functions

        // Access

            //- Reference to mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }
            fvMesh& mesh()
            {
                return mesh_;
            }

            scalar mergeDistance() const
            {
                return mergeDistance_;
            }

            //- Overwrite the mesh?
            bool overwrite() const
            {
                return overwrite_;
            }

            //- (points)instance of mesh upon construction
            const word& oldInstance() const
            {
                return oldInstance_;
            }

            //- Reference to surface search engines
            const refinementSurfaces& surfaces() const
            {
                return surfaces_;
            }

            //- Reference to feature edge mesh
            const refinementFeatures& features() const
            {
                return features_;
            }

            //- Reference to refinement shells (regions)
            const refinementRegions& shells() const
            {
                return shells_;
            }

            //- Reference to meshcutting engine
            const hexRef8& meshCutter() const
            {
                return meshCutter_;
            }

            //- Per start-end edge the index of the surface hit
            const labelList& surfaceIndex() const
            {
                return surfaceIndex_;
            }

            labelList& surfaceIndex()
            {
                return surfaceIndex_;
            }

            //- Additional face data that is maintained across
            //  topo changes. Every entry is a list over all faces.
            //  Bit of a hack. Additional flag to say whether to maintain master
            //  only (false) or increase set to account for face-from-face.
            const List<Tuple2<mapType, labelList>>& userFaceData() const
            {
                return userFaceData_;
            }

            List<Tuple2<mapType, labelList>>& userFaceData()
            {
                return userFaceData_;
            }


        // Other

            //- Count number of intersections (local)
            label countHits() const;

            //- Redecompose according to cell count
            //  keepZoneFaces : find all faceZones from zoned surfaces and keep
            //                  owner and neighbour together
            //  keepBaffles   : find all baffles and keep them together
            autoPtr<mapDistributePolyMesh> balance
            (
                const bool keepZoneFaces,
                const bool keepBaffles,
                const scalarField& cellWeights,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor
            );

            //- Get faces with intersection.
            labelList intersectedFaces() const;

            //- Get points on surfaces with intersection and boundary faces.
            labelList intersectedPoints() const;

            //- Create patch from set of patches
            static autoPtr<indirectPrimitivePatch> makePatch
            (
                const polyMesh&,
                const labelList&
            );

            //- Helper function to make a pointVectorField with correct
            //  bcs for mesh movement:
            //  - adaptPatchIDs         : fixedValue
            //  - processor             : calculated (so free to move)
            //  - cyclic/wedge/symmetry : slip
            //  - other                 : slip
            static tmp<pointVectorField> makeDisplacementField
            (
                const pointMesh& pMesh,
                const labelList& adaptPatchIDs
            );

            //- Helper function: check that face zones are synced
            static void checkCoupledFaceZones(const polyMesh&);

            //- Helper: calculate edge weights (1/length)
            static void calculateEdgeWeights
            (
                const polyMesh& mesh,
                const PackedBoolList& isMasterEdge,
                const labelList& meshPoints,
                const edgeList& edges,
                scalarField& edgeWeights,
                scalarField& invSumWeight
            );

            //- Helper: weighted sum (over all subset of mesh points) by
            //  summing contribution from (master) edges
            template<class Type>
            static void weightedSum
            (
                const polyMesh& mesh,
                const PackedBoolList& isMasterEdge,
                const labelList& meshPoints,
                const edgeList& edges,
                const scalarField& edgeWeights,
                const Field<Type>& data,
                Field<Type>& sum
            );


        // Refinement

            //- Is local topology a small gap?
            bool isGap
            (
                const scalar,
                const vector&,
                const vector&,
                const vector&,
                const vector&
            ) const;

            //- Is local topology a small gap normal to the test vector
            bool isNormalGap
            (
                const scalar,
                const vector&,
                const vector&,
                const vector&,
                const vector&
            ) const;

            //- Calculate list of cells to refine.
            labelList refineCandidates
            (
                const pointField& keepPoints,
                const scalar curvature,
                const scalar planarAngle,

                const bool featureRefinement,
                const bool featureDistanceRefinement,
                const bool internalRefinement,
                const bool surfaceRefinement,
                const bool curvatureRefinement,
                const bool gapRefinement,
                const label maxGlobalCells,
                const label maxLocalCells
            ) const;

            //- Refine some cells
            autoPtr<mapPolyMesh> refine(const labelList& cellsToRefine);

            //- Refine some cells and rebalance
            autoPtr<mapDistributePolyMesh> refineAndBalance
            (
                const string& msg,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor,
                const labelList& cellsToRefine,
                const scalar maxLoadUnbalance
            );

            //- Balance before refining some cells
            autoPtr<mapDistributePolyMesh> balanceAndRefine
            (
                const string& msg,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor,
                const labelList& cellsToRefine,
                const scalar maxLoadUnbalance
            );


        // Baffle handling

            //- Split off unreachable areas of mesh.
            void baffleAndSplitMesh
            (
                const bool handleSnapProblems,

                // How to remove problem snaps
                const snapParameters& snapParams,
                const bool useTopologicalSnapDetection,
                const bool removeEdgeConnectedCells,
                const scalarField& perpendicularAngle,

                // How to handle free-standing baffles
                const bool mergeFreeStanding,
                const scalar freeStandingAngle,

                const dictionary& motionDict,
                Time& runTime,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,
                const point& keepPoint
            );

            //- Split off (with optional buffer layers) unreachable areas
            //  of mesh. Does not introduce baffles.
            autoPtr<mapPolyMesh> splitMesh
            (
                const label nBufferLayers,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,
                const point& keepPoint
            );

            //- Find boundary points that connect to more than one cell
            //  region and split them.
            autoPtr<mapPolyMesh> dupNonManifoldPoints(const localPointRegion&);

            //- Find boundary points that connect to more than one cell
            //  region and split them.
            autoPtr<mapPolyMesh> dupNonManifoldPoints();

            //- Create baffle for every internal face where ownPatch != -1.
            //  External faces get repatched according to ownPatch (nbrPatch
            //  should be -1 for these)
            autoPtr<mapPolyMesh> createBaffles
            (
                const labelList& ownPatch,
                const labelList& nbrPatch
            );

            //- Debug helper: check faceZones are not on processor patches
            void checkZoneFaces() const;

            //- Create baffles for faces straddling zoned surfaces. Return
            //  baffles.
            autoPtr<mapPolyMesh> createZoneBaffles
            (
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,
                List<labelPair>&
            );

            //- Merge baffles. Gets pairs of faces.
            autoPtr<mapPolyMesh> mergeBaffles(const List<labelPair>&);

            //- Put faces/cells into zones according to surface specification.
            //  Returns null if no zone surfaces present. Region containing
            //  the keepPoint will not be put into a cellZone.
            autoPtr<mapPolyMesh> zonify
            (
                const point& keepPoint,
                const bool allowFreeStandingZoneFaces
            );


        // Other topo changes

            //- Add patch originating from meshing. Update meshedPatches_.
            label addMeshedPatch(const word& name, const dictionary&);

            //- Get patchIDs for patches added in addMeshedPatch.
            labelList meshedPatches() const;

            //- Select coupled faces that are not collocated
            void selectSeparatedCoupledFaces(boolList&) const;

            //- Find region point is in. Uses optional perturbation to re-test.
            static label findRegion
            (
                const polyMesh&,
                const labelList& cellRegion,
                const vector& perturbVec,
                const point& p
            );

            //- Split mesh. Keep part containing point.
            autoPtr<mapPolyMesh> splitMeshRegions
            (
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,
                const point& keepPoint
            );

            //- Split faces into two
            autoPtr<mapPolyMesh> splitFaces
            (
                const labelList& splitFaces,
                const labelPairList& splits
            );

            //- Update local numbering for mesh redistribution
            void distribute(const mapDistributePolyMesh&);

            //- Update for external change to mesh. changedFaces are in new mesh
            //  face labels.
            void updateMesh
            (
                const mapPolyMesh&,
                const labelList& changedFaces
            );

            //- Helper: reorder list according to map.
            template<class T>
            static void updateList
            (
                const labelList& newToOld,
                const T& nullValue,
                List<T>& elems
            );


            // Restoring : is where other processes delete and reinsert data.

                //- Signal points/face/cells for which to store data
                void storeData
                (
                    const labelList& pointsToStore,
                    const labelList& facesToStore,
                    const labelList& cellsToStore
                );

                //- Update local numbering + undo
                //  Data to restore given as new pointlabel + stored pointlabel
                //  (i.e. what was in pointsToStore)
                void updateMesh
                (
                    const mapPolyMesh&,
                    const labelList& changedFaces,
                    const Map<label>& pointsToRestore,
                    const Map<label>& facesToRestore,
                    const Map<label>& cellsToRestore
                );

            // Merging coplanar faces and edges

                //- Merge coplanar faces. preserveFaces is != -1 for faces
                //  to be preserved
                label mergePatchFacesUndo
                (
                    const scalar minCos,
                    const scalar concaveCos,
                    const labelHashSet& patchIDs,
                    const dictionary& motionDict,
                    const labelList& preserveFaces
                );

                autoPtr<mapPolyMesh> doRemovePoints
                (
                    removePoints& pointRemover,
                    const boolList& pointCanBeDeleted
                );

                autoPtr<mapPolyMesh> doRestorePoints
                (
                    removePoints& pointRemover,
                    const labelList& facesToRestore
                );

                labelList collectFaces
                (
                    const labelList& candidateFaces,
                    const labelHashSet& set
                ) const;

                // Pick up faces of cells of faces in set.
                labelList growFaceCellFace
                (
                    const labelHashSet& set
                ) const;

                //- Merge edges, maintain mesh quality. Return global number
                //  of edges merged
                label mergeEdgesUndo
                (
                    const scalar minCos,
                    const dictionary& motionDict
                );


        // Debug/IO

            //- Debugging: check that all faces still obey start()>end()
            void checkData();

            static void testSyncPointList
            (
                const string& msg,
                const polyMesh& mesh,
                const List<scalar>& fld
            );

            static void testSyncPointList
            (
                const string& msg,
                const polyMesh& mesh,
                const List<point>& fld
            );

            //- Compare two lists over all boundary faces
            template<class T>
            void testSyncBoundaryFaceList
            (
                const scalar mergeDistance,
                const string&,
                const UList<T>&,
                const UList<T>&
            ) const;

            //- Print list according to (collected and) sorted coordinate
            template<class T>
            static void collectAndPrint
            (
                const UList<point>& points,
                const UList<T>& data
            );

            //- Determine master point for subset of points. If coupled
            //  chooses only one
            static PackedBoolList getMasterPoints
            (
                const polyMesh& mesh,
                const labelList& meshPoints
            );

            //- Determine master edge for subset of edges. If coupled
            //  chooses only one
            static PackedBoolList getMasterEdges
            (
                const polyMesh& mesh,
                const labelList& meshEdges
            );

            //- Print some mesh stats.
            void printMeshInfo(const bool, const string&) const;

            //- Replacement for Time::timeName() : return oldInstance (if
            //  overwrite_)
            word timeName() const;

            //- Set instance of all local IOobjects
            void setInstance(const fileName&);

            //- Write mesh and all data
            bool write() const;

            //- Write refinement level as volScalarFields for postprocessing
            void dumpRefinementLevel() const;

            //- Debug: Write intersection information to OBJ format
            void dumpIntersections(const fileName& prefix) const;

            //- Do any one of above IO functions
            void write
            (
                const debugType debugFlags,
                const writeType writeFlags,
                const fileName&
            ) const;

            //- Helper: calculate average
            template<class T>
            static T gAverage
            (
                const PackedBoolList& isMasterElem,
                const UList<T>& values
            );

            //- Get/set write level
            static writeType writeLevel();
            static void writeLevel(const writeType);

            //- Get/set output level
            static outputType outputLevel();
            static void outputLevel(const outputType);


            //- Helper: convert wordList into bit pattern using provided
            //  NamedEnum
            template<class Enum>
            static int readFlags(const Enum& namedEnum, const wordList&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshRefinement&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshRefinementTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
