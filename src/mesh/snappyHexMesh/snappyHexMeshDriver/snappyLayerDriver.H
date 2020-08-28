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
    Foam::snappyLayerDriver

Description
    All to do with adding layers

SourceFiles
    snappyLayerDriver.C

\*---------------------------------------------------------------------------*/

#ifndef snappyLayerDriver_H
#define snappyLayerDriver_H

#include "meshRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class removePoints;
class pointSet;
class motionSmoother;
class addPatchCellLayer;
class pointData;
class faceSet;
class layerParameters;

/*---------------------------------------------------------------------------*\
                      Class snappyLayerDriver Declaration
\*---------------------------------------------------------------------------*/

class snappyLayerDriver
{
public:

    // Public data types

        //- Extrusion controls
        enum extrudeMode
        {
            NOEXTRUDE,      //!< Do not extrude. No layers added.
            EXTRUDE,        //!< Extrude
            EXTRUDEREMOVE   //!< Extrude but afterwards remove added
                            //!< faces locally
        };

private:

    // Private classes

        //- Combine operator class to combine normal with other normal.
        class normalsCombine
        {
        public:

            void operator()(vector& x, const vector& y) const
            {
                if (y != point::max)
                {
                    if (x == point::max)
                    {
                        x = y;
                    }
                    else
                    {
                        x *= (x&y);
                    }
                }
            }
        };


    // Private Data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- From surface region to patch
        const labelList globalToMasterPatch_;

        //- From surface region to patch
        const labelList globalToSlavePatch_;



    // Private Member Functions

        // Layers

            //- For debugging: Dump displacement to .obj files
            static void dumpDisplacement
            (
                const fileName&,
                const indirectPrimitivePatch&,
                const vectorField&,
                const List<extrudeMode>&
            );

            //- Average point wise data to face wise
            static tmp<scalarField> avgPointData
            (
                const indirectPrimitivePatch&,
                const scalarField& pointFld
            );

            //- Check that primitivePatch is not multiply connected.
            //  Collect non-manifold points in pointSet.
            static void checkManifold
            (
                const indirectPrimitivePatch&,
                pointSet& nonManifoldPoints
            );

            //- Check that mesh outside is not multiply connected.
            void checkMeshManifold() const;


            // Static extrusion setup

                //- Unset extrusion on point. Returns true if anything unset.
                static bool unmarkExtrusion
                (
                    const label patchPointi,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- Unset extrusion on face. Returns true if anything unset.
                static bool unmarkExtrusion
                (
                    const face& localFace,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- No extrusion at non-manifold points.
                void handleNonManifolds
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,
                    const labelListList& edgeGlobalFaces,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- No extrusion on feature edges. Assumes non-manifold
                //  edges already handled.
                void handleFeatureAngle
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,
                    const scalar minCos,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- No extrusion on warped faces
                void handleWarpedFaces
                (
                    const indirectPrimitivePatch& pp,
                    const scalar faceRatio,
                    const scalar edge0Len,
                    const labelList& cellLevel,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Determine the number of layers per point from the number of
                //  layers per surface.
                void setNumLayers
                (
                    const labelList& patchToNLayers,
                    const labelList& patchIDs,
                    const indirectPrimitivePatch& pp,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus,
                    label& nIdealAddedCells
                ) const;

                //- Helper function to make a pointVectorField with correct
                //  bcs for layer addition:
                //  - numLayers > 0         : fixedValue
                //  - numLayers == 0        : fixedValue (always zero)
                //  - processor             : calculated (so free to move)
                //  - cyclic/wedge/symmetry : slip
                //  - other                 : slip
                static tmp<pointVectorField> makeLayerDisplacementField
                (
                    const pointMesh& pMesh,
                    const labelList& numLayers
                );

                //- Grow no-extrusion layer.
                void growNoExtrusion
                (
                    const indirectPrimitivePatch& pp,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- See what patches boundaryedges should be extruded into
                void determineSidePatches
                (
                    const globalIndex& globalFaces,
                    const labelListList& edgeGlobalFaces,
                    const indirectPrimitivePatch& pp,

                    labelList& sidePatchID
                );

                //- Calculate pointwise wanted and minimum thickness.
                //  thickness: wanted thickness
                //  minthickness: when to give up and not extrude
                //  Gets per patch parameters and determine pp pointwise
                //  parameters.
                void calculateLayerThickness
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& patchIDs,
                    const layerParameters& layerParams,
                    const labelList& cellLevel,
                    const labelList& patchNLayers,
                    const scalar edge0Len,

                    scalarField& thickness,
                    scalarField& minThickness,
                    scalarField& expansionRatio
                ) const;


            // Extrusion execution

                //- Synchronize displacement among coupled patches.
                void syncPatchDisplacement
                (
                    const indirectPrimitivePatch& pp,
                    const scalarField& minThickness,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Get nearest point on surface to snap to
                void getPatchDisplacement
                (
                    const indirectPrimitivePatch& pp,
                    const scalarField& thickness,
                    const scalarField& minThickness,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- For truncateDisplacement: find strings of edges
                bool sameEdgeNeighbour
                (
                    const labelListList& globalEdgeFaces,
                    const label myGlobalFacei,
                    const label nbrGlobFacei,
                    const label edgeI
                ) const;

                //- For truncateDisplacement: find strings of edges
                void getVertexString
                (
                    const indirectPrimitivePatch& pp,
                    const labelListList& globalEdgeFaces,
                    const label facei,
                    const label edgeI,
                    const label myGlobFacei,
                    const label nbrGlobFacei,
                    DynamicList<label>& vertices
                ) const;

                //- Truncates displacement
                // - for all patchFaces in the faceset displacement gets set
                //   to zero
                // - all displacement < minThickness gets set to zero
                // - all non-consecutive extrusions get set to 0
                label truncateDisplacement
                (
                    const globalIndex& globalFaces,
                    const labelListList& edgeGlobalFaces,
                    const indirectPrimitivePatch& pp,
                    const scalarField& minThickness,
                    const faceSet& illegalPatchFaces,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Setup layer information (at points and faces) to
                //  modify mesh topology in
                //  regions where layer mesh terminates. Guarantees an
                //  optional slow decreasing of the number of layers.
                //  Returns the number of layers per face and per point
                //  to go into the actual layer addition engine.
                void setupLayerInfoTruncation
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& patchNLayers,
                    const List<extrudeMode>& extrudeStatus,
                    const label nBufferCellsNoExtrude,
                    labelList& nPatchPointLayers,
                    labelList& nPatchFaceLayers
                ) const;

                //- Does any of the cells use a face from faces?
                static bool cellsUseFace
                (
                    const polyMesh& mesh,
                    const labelList& cellLabels,
                    const labelHashSet& faces
                );

                //- Checks the newly added cells and locally unmarks points
                //  so they will not get extruded next time round. Returns
                //  global number of unmarked points (0 if all was fine)
                static label checkAndUnmark
                (
                    const addPatchCellLayer& addLayer,
                    const dictionary& motionDict,
                    const bool additionalReporting,
                    const List<labelPair>& baffles,
                    const indirectPrimitivePatch& pp,
                    const fvMesh&,

                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- Count global number of extruded faces
                static label countExtrusion
                (
                    const indirectPrimitivePatch& pp,
                    const List<extrudeMode>& extrudeStatus
                );

                //- Collect layer faces and layer cells into bools
                //  for ease of handling
                static void getLayerCellsFaces
                (
                    const polyMesh&,
                    const addPatchCellLayer&,
                    const scalarField& oldRealThickness,

                    labelList& cellStatus,
                    scalarField& faceRealThickness
                );

                //- Print layer coverage table
                void printLayerData
                (
                    const fvMesh& mesh,
                    const labelList& patchIDs,
                    const labelList& cellNLayers,
                    const scalarField& faceWantedThickness,
                    const scalarField& faceRealThickness
                ) const;

                //- Write cellSet,faceSet for layers
                bool writeLayerData
                (
                    const fvMesh& mesh,
                    const labelList& patchIDs,
                    const labelList& cellNLayers,
                    const scalarField& faceWantedThickness,
                    const scalarField& faceRealThickness
                ) const;


            // Mesh shrinking (to create space for layers)

                //- Average field (over all subset of mesh points) by
                //  summing contribution from edges. Global parallel since only
                //  does master edges for coupled edges.
                template<class Type>
                static void averageNeighbours
                (
                    const polyMesh& mesh,
                    const PackedBoolList& isMasterEdge,
                    const labelList& meshEdges,
                    const labelList& meshPoints,
                    const edgeList& edges,
                    const scalarField& invSumWeight,
                    const Field<Type>& data,
                    Field<Type>& average
                );

                //- Calculate inverse sum of edge weights (currently always 1.0)
                void sumWeights
                (
                    const PackedBoolList& isMasterEdge,
                    const labelList& meshEdges,
                    const labelList& meshPoints,
                    const edgeList& edges,
                    scalarField& invSumWeight
                ) const;

                //- Smooth scalar field on patch
                void smoothField
                (
                    const motionSmoother& meshMover,
                    const PackedBoolList& isMasterPoint,
                    const PackedBoolList& isMasterEdge,
                    const labelList& meshEdges,
                    const scalarField& fieldMin,
                    const label nSmoothDisp,
                    scalarField& field
                ) const;

                //- Smooth normals on patch.
                void smoothPatchNormals
                (
                    const motionSmoother& meshMover,
                    const PackedBoolList& isMasterPoint,
                    const PackedBoolList& isMasterEdge,
                    const labelList& meshEdges,
                    const label nSmoothDisp,
                    pointField& normals
                ) const;

                //- Smooth normals in interior.
                void smoothNormals
                (
                    const label nSmoothDisp,
                    const PackedBoolList& isMasterPoint,
                    const PackedBoolList& isMasterEdge,
                    const labelList& fixedPoints,
                    pointVectorField& normals
                ) const;

                bool isMaxEdge
                (
                    const List<pointData>&,
                    const label edgeI,
                    const scalar minCos
                ) const;

                //- Stop layer growth where mesh wraps around edge with a
                //  large feature angle
                void handleFeatureAngleLayerTerminations
                (
                    const scalar minCos,
                    const PackedBoolList& isMasterPoint,
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,

                    List<extrudeMode>& extrudeStatus,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    label& nPointCounter
                ) const;

                //- Find isolated islands (points, edges and faces and
                // layer terminations)
                // in the layer mesh and stop any layer growth at these points.
                void findIsolatedRegions
                (
                    const scalar minCosLayerTermination,
                    const PackedBoolList& isMasterPoint,
                    const PackedBoolList& isMasterEdge,
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,
                    const scalarField& minThickness,

                    List<extrudeMode>& extrudeStatus,
                    pointField& patchDisp,
                    labelList& patchNLayers
                ) const;

                // Calculate medial axis fields
                void medialAxisSmoothingInfo
                (
                    const motionSmoother& meshMover,
                    const label nSmoothNormals,
                    const label nSmoothSurfaceNormals,
                    const scalar minMedialAxisAngleCos,
                    const scalar featureAngle,

                    pointVectorField& dispVec,
                    pointScalarField& medialRatio,
                    pointScalarField& medialDist,
                    pointVectorField& medialVec
                ) const;

                //- Main routine to shrink mesh
                void shrinkMeshMedialDistance
                (
                    motionSmoother& meshMover,
                    const dictionary& meshQualityDict,
                    const List<labelPair>& baffles,
                    const label nSmoothPatchThickness,
                    const label nSmoothDisplacement,
                    const scalar maxThicknessToMedialRatio,
                    const label nAllowableErrors,
                    const label nSnap,
                    const scalar minCosLayerTermination,

                    const scalarField& layerThickness,
                    const scalarField& minThickness,

                    const pointVectorField& dispVec,
                    const pointScalarField& medialRatio,
                    const pointScalarField& medialDist,
                    const pointVectorField& medialVec,

                    List<extrudeMode>& extrudeStatus,
                    pointField& patchDisp,
                    labelList& patchNLayers
                ) const;


public:

    //- Runtime type information
    ClassName("snappyLayerDriver");

    // Constructors

        //- Construct from components
        snappyLayerDriver
        (
            meshRefinement& meshRefiner,
            const labelList& globalToMasterPatch,
            const labelList& globalToSlavePatch
        );

        //- Disallow default bitwise copy construction
        snappyLayerDriver(const snappyLayerDriver&) = delete;


    // Member Functions

            //- Merge patch faces on same cell.
            void mergePatchFacesUndo
            (
                const layerParameters& layerParams,
                const dictionary& motionDict
            );

            //- Add cell layers
            void addLayers
            (
                const layerParameters& layerParams,
                const dictionary& motionDict,
                const labelList& patchIDs,
                const label nAllowableErrors,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor
            );

            //- Add layers according to the dictionary settings
            void doLayers
            (
                const dictionary& shrinkDict,
                const dictionary& motionDict,
                const layerParameters& layerParams,
                const bool preBalance,              // balance before adding?
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const snappyLayerDriver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "snappyLayerDriverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
