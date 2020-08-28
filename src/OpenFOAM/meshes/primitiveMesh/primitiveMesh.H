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
    Foam::primitiveMesh

Description
    Cell-face mesh analysis engine

SourceFiles
    primitiveMeshI.H
    primitiveMesh.C
    primitiveMeshClear.C
    primitiveMeshCellCells.C
    primitiveMeshEdgeCells.C
    primitiveMeshPointCells.C
    primitiveMeshCells.C
    primitiveMeshEdgeFaces.C
    primitiveMeshPointFaces.C
    primitiveMeshCellEdges.C
    primitiveMeshPointEdges.C
    primitiveMeshPointPoints.C
    primitiveMeshEdges.C
    primitiveMeshCellCentresAndVols.C
    primitiveMeshFaceCentresAndAreas.C
    primitiveMeshFindCell.C

\*---------------------------------------------------------------------------*/

#ifndef primitiveMesh_H
#define primitiveMesh_H

#include "DynamicList.H"
#include "edgeList.H"
#include "pointField.H"
#include "faceList.H"
#include "cellList.H"
#include "cellShapeList.H"
#include "labelList.H"
#include "boolList.H"
#include "HashSet.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class PackedBoolList;

/*---------------------------------------------------------------------------*\
                        Class primitiveMesh Declaration
\*---------------------------------------------------------------------------*/

class primitiveMesh
{
    // Permanent data

        // Primitive size data

            //- Number of internal points (or -1 if points not sorted)
            label nInternalPoints_;

            //- Number of points
            label nPoints_;

            //- Number of internal edges using 0 boundary points
            mutable label nInternal0Edges_;

            //- Number of internal edges using 0 or 1 boundary points
            mutable label nInternal1Edges_;

            //- Number of internal edges using 0,1 or 2boundary points
            mutable label nInternalEdges_;

            //- Number of edges
            mutable label nEdges_;

            //- Number of internal faces
            label nInternalFaces_;

            //- Number of faces
            label nFaces_;

            //- Number of cells
            label nCells_;


        // Shapes

            //- Cell shapes
            mutable cellShapeList* cellShapesPtr_;

            //- Edges
            mutable edgeList* edgesPtr_;


        // Connectivity

            //- Cell-cells
            mutable labelListList* ccPtr_;

            //- Edge-cells
            mutable labelListList* ecPtr_;

            //- Point-cells
            mutable labelListList* pcPtr_;

            //- Cell-faces
            mutable cellList* cfPtr_;

            //- Edge-faces
            mutable labelListList* efPtr_;

            //- Point-faces
            mutable labelListList* pfPtr_;

            //- Cell-edges
            mutable labelListList* cePtr_;

            //- Face-edges
            mutable labelListList* fePtr_;

            //- Point-edges
            mutable labelListList* pePtr_;

            //- Point-points
            mutable labelListList* ppPtr_;

            //- Cell-points
            mutable labelListList* cpPtr_;


        // On-the-fly edge addressing storage

            //- Temporary storage for addressing.
            mutable DynamicList<label> labels_;

            //- Temporary storage for addressing
            mutable labelHashSet labelSet_;


        // Geometric data

            //- Cell centres
            mutable vectorField* cellCentresPtr_;

            //- Face centres
            mutable vectorField* faceCentresPtr_;

            //- Cell volumes
            mutable scalarField* cellVolumesPtr_;

            //- Face areas
            mutable vectorField* faceAreasPtr_;

            //- Face area magnitudes
            mutable scalarField* magFaceAreasPtr_;


        // Topological calculations

            //- Calculate cell shapes
            void calcCellShapes() const;

            //- Calculate cell-cell addressing
            void calcCellCells() const;

            //- Calculate point-cell addressing
            void calcPointCells() const;

            //- Calculate cell-face addressing
            void calcCells() const;

            //- Calculate edge list
            void calcCellEdges() const;

            //- Calculate point-point addressing
            void calcPointPoints() const;

            //- Calculate edges, pointEdges and faceEdges (if doFaceEdges=true)
            //  During edge calculation, a larger set of data is assembled.
            //  Create and destroy as a set, using clearOutEdges()
            void calcEdges(const bool doFaceEdges) const;
            void clearOutEdges();
            //- Helper: return (after optional creation) edge between two points
            static label getEdge
            (
                List<DynamicList<label>>&,
                DynamicList<edge>&,
                const label,
                const label
            );
            //- For on-the-fly addressing calculation
            static label findFirstCommonElementFromSortedLists
            (
                const labelList&,
                const labelList&
            );

protected:

    // Static Data Members

        //- Static data to control mesh checking

            //- Cell closedness warning threshold
            //  set as the fraction of un-closed area to closed area
            static scalar closedThreshold_;

            //- Aspect ratio warning threshold
            static scalar aspectThreshold_;

            //- Non-orthogonality warning threshold in deg
            static scalar nonOrthThreshold_;

            //- Skewness warning threshold
            static scalar skewThreshold_;

            //- Threshold where faces are considered coplanar
            static scalar planarCosAngle_;


        // Geometrical calculations

            //- Calculate face centres and areas
            void calcFaceCentresAndAreas() const;

            void makeFaceCentresAndAreas
            (
                const pointField& p,
                vectorField& fCtrs,
                vectorField& fAreas,
                scalarField& magfAreas
            ) const;

            //- Calculate cell centres and volumes
            void calcCellCentresAndVols() const;

            void makeCellCentresAndVols
            (
                const vectorField& fCtrs,
                const vectorField& fAreas,
                vectorField& cellCtrs,
                scalarField& cellVols
            ) const;

            //- Calculate edge vectors
            void calcEdgeVectors() const;


        // Mesh checking

            //- Check if all points on face are shared with another face.
            bool checkDuplicateFaces
            (
                const label,
                const Map<label>&,
                label& nBaffleFaces,
                labelHashSet*
            ) const;

            //- Check that shared points are in consecutive order.
            bool checkCommonOrder
            (
                const label,
                const Map<label>&,
                labelHashSet*
            ) const;

            //- Check boundary for closedness
            bool checkClosedBoundary
            (
                const vectorField&,
                const bool,
                const PackedBoolList&
            ) const;

            //- Check cells for closedness
            bool checkClosedCells
            (
                const vectorField& faceAreas,
                const scalarField& cellVolumes,
                const bool report,
                labelHashSet* setPtr,
                labelHashSet* aspectSetPtr,
                const Vector<label>& meshD
            ) const;

            //- Check for negative face areas
            bool checkFaceAreas
            (
                const vectorField& faceAreas,
                const bool report,
                const bool detailedReport,
                labelHashSet* setPtr
            ) const;

            //- Check for negative cell volumes
            bool checkCellVolumes
            (
                const scalarField& vols,
                const bool report,
                const bool detailedReport,
                labelHashSet* setPtr
            ) const;

            //- Check for non-orthogonality
            bool checkFaceOrthogonality
            (
                const vectorField& fAreas,
                const vectorField& cellCtrs,
                const bool report,
                labelHashSet* setPtr
            ) const;

            //- Check face pyramid volume
            bool checkFacePyramids
            (
                const pointField& points,
                const vectorField& ctrs,
                const bool report,
                const bool detailedReport,
                const scalar minPyrVol,
                labelHashSet* setPtr
            ) const;

            //- Check face skewness
            bool checkFaceSkewness
            (
                const pointField& points,
                const vectorField& fCtrs,
                const vectorField& fAreas,
                const vectorField& cellCtrs,
                const bool report,
                labelHashSet* setPtr
            ) const;

            //- Check face angles
            //  Allows a slight non-convexity.  E.g. maxDeg = 10 allows for
            //  angles < 190 (or 10 degrees concavity) (if truly concave and
            //  points not visible from face centre the face-pyramid check in
            //  checkMesh will fail)
            bool checkFaceAngles
            (
                const pointField& points,
                const vectorField& faceAreas,
                const bool report,
                const scalar maxDeg,
                labelHashSet* setPtr
            ) const;

            //- Check face warpage
            bool checkFaceFlatness
            (
                const pointField& points,
                const vectorField& faceCentres,
                const vectorField& faceAreas,
                const bool report,
                const scalar warnFlatness,
                labelHashSet* setPtr
            ) const;

            //- Check for concave cells by the planes of faces
            bool checkConcaveCells
            (
                const vectorField& fAreas,
                const pointField& fCentres,
                const bool report,
                labelHashSet* setPtr
            ) const;


        //- Construct null
        primitiveMesh();


public:

        // Static data

            ClassName("primitiveMesh");

            //- Estimated number of cells per edge
            static const unsigned cellsPerEdge_ = 4;

            //- Estimated number of cells per point
            static const unsigned cellsPerPoint_ = 8;

            //- Estimated number of faces per cell
            static const unsigned facesPerCell_ = 6;

            //- Estimated number of faces per edge
            static const unsigned facesPerEdge_ = 4;

            //- Estimated number of faces per point
            static const unsigned facesPerPoint_ = 12;

            //- Estimated number of edges per cell
            static const unsigned edgesPerCell_ = 12;

            //- Estimated number of edges per cell
            static const unsigned edgesPerFace_ = 4;

            //- Estimated number of edges per point
            static const unsigned edgesPerPoint_ = 6;

            //- Estimated number of points per cell
            static const unsigned pointsPerCell_ = 8;

            //- Estimated number of points per face
            static const unsigned pointsPerFace_ = 4;


    // Constructors

        //- Construct from components
        primitiveMesh
        (
            const label nPoints,
            const label nInternalFaces,
            const label nFaces,
            const label nCells
        );

        //- Disallow default bitwise copy construction
        primitiveMesh(const primitiveMesh&);


    //- Destructor
    virtual ~primitiveMesh();


    // Member Functions

        //- Reset this primitiveMesh given the primitive array sizes
        void reset
        (
            const label nPoints,
            const label nInternalFaces,
            const label nFaces,
            const label nCells
        );

        //- Reset this primitiveMesh given the primitive array sizes and cells
        void reset
        (
            const label nPoints,
            const label nInternalFaces,
            const label nFaces,
            const label nCells,
            cellList& cells
        );

        //- Reset this primitiveMesh given the primitive array sizes and cells
        void reset
        (
            const label nPoints,
            const label nInternalFaces,
            const label nFaces,
            const label nCells,
            cellList&& cells
        );


        // Access

            // Mesh size parameters

                inline label nPoints() const;
                inline label nEdges() const;
                inline label nInternalFaces() const;
                inline label nFaces() const;
                inline label nCells() const;

                // If points are ordered (nInternalPoints != -1):

                    //- Points not on boundary
                    inline label nInternalPoints() const;

                    //- Internal edges (i.e. not on boundary face) using
                    //  no boundary point
                    inline label nInternal0Edges() const;

                    //- Internal edges using 0 or 1 boundary point
                    inline label nInternal1Edges() const;

                    //- Internal edges using 0,1 or 2 boundary points
                    inline label nInternalEdges() const;


            // Primitive mesh data

                //- Return mesh points
                virtual const pointField& points() const = 0;

                //- Return faces
                virtual const faceList& faces() const = 0;

                //- Face face-owner addressing
                virtual const labelList& faceOwner() const = 0;

                //- Face face-neighbour addressing
                virtual const labelList& faceNeighbour() const = 0;

                //- Return old points for mesh motion
                virtual const pointField& oldPoints() const = 0;


            // Derived mesh data

                //- Return cell shapes
                const cellShapeList& cellShapes() const;

                //- Return mesh edges. Uses calcEdges.
                const edgeList& edges() const;

                //- Helper function to calculate cell-face addressing from
                //  face-cell addressing. If nCells is not provided it will
                //  scan for the maximum.
                static void calcCells
                (
                    cellList&,
                    const labelUList& own,
                    const labelUList& nei,
                    const label nCells = -1
                );

                //- Helper function to calculate point ordering. Returns true
                //  if points already ordered, false and fills pointMap (old to
                //  new). Map splits points into those not used by any boundary
                //  face and those that are.
                static bool calcPointOrder
                (
                    label& nInternalPoints,
                    labelList& pointMap,
                    const faceList&,
                    const label nInternalFaces,
                    const label nPoints
                );

            // Return mesh connectivity

                const labelListList& cellCells() const;
                // faceCells given as owner and neighbour
                const labelListList& edgeCells() const;
                const labelListList& pointCells() const;

                const cellList& cells() const;
                // faceFaces considered unnecessary
                const labelListList& edgeFaces() const;
                const labelListList& pointFaces() const;

                const labelListList& cellEdges() const;
                const labelListList& faceEdges() const;
                // edgeEdges considered unnecessary
                const labelListList& pointEdges() const;
                const labelListList& pointPoints() const;
                const labelListList& cellPoints() const;


            // Geometric data (raw!)

                const vectorField& cellCentres() const;
                const vectorField& faceCentres() const;
                const scalarField& cellVolumes() const;
                const vectorField& faceAreas() const;
                const scalarField& magFaceAreas() const;


            // Mesh motion

                //- Move points, returns volumes swept by faces in motion
                tmp<scalarField> movePoints
                (
                    const pointField& p,
                    const pointField& oldP
                );


            //- Return true if given face label is internal to the mesh
            inline bool isInternalFace(const label faceIndex) const;


            // Topological checks

                //- Check face ordering
                virtual bool checkUpperTriangular
                (
                    const bool report = false,
                    labelHashSet* setPtr = nullptr
                ) const;

                //- Check cell zip-up
                virtual bool checkCellsZipUp
                (
                    const bool report = false,
                    labelHashSet* setPtr = nullptr
                ) const;

                //- Check uniqueness of face vertices
                virtual bool checkFaceVertices
                (
                    const bool report = false,
                    labelHashSet* setPtr = nullptr
                ) const;

                //- Check for unused points
                virtual bool checkPoints
                (
                    const bool report = false,
                    labelHashSet* setPtr = nullptr
                ) const;

                //- Check face-face connectivity
                virtual bool checkFaceFaces
                (
                    const bool report = false,
                    labelHashSet* setPtr = nullptr
                ) const;


            // Geometric checks

                //- Check boundary for closedness
                virtual bool checkClosedBoundary(const bool report = false)
                const;

                //- Check cells for closedness
                virtual bool checkClosedCells
                (
                    const bool report = false,
                    labelHashSet* setPtr = nullptr,
                    labelHashSet* highAspectSetPtr = nullptr,
                    const Vector<label>& solutionD = Vector<label>::one
                ) const;

                //- Check for negative face areas
                virtual bool checkFaceAreas
                (
                    const bool report = false,
                    labelHashSet* setPtr = nullptr
                ) const;

                //- Check for negative cell volumes
                virtual bool checkCellVolumes
                (
                    const bool report = false,
                    labelHashSet* setPtr = nullptr
                ) const;

                //- Check for non-orthogonality
                virtual bool checkFaceOrthogonality
                (
                    const bool report = false,
                    labelHashSet* setPtr = nullptr
                ) const;

                //- Check face pyramid volume
                virtual bool checkFacePyramids
                (
                    const bool report = false,
                    const scalar minPyrVol = -small,
                    labelHashSet* setPtr = nullptr
                ) const;

                //- Check face skewness
                virtual bool checkFaceSkewness
                (
                    const bool report = false,
                    labelHashSet* setPtr = nullptr
                ) const;

                //- Check face angles
                virtual bool checkFaceAngles
                (
                    const bool report = false,
                    const scalar maxSin = 10,    // In degrees
                    labelHashSet* setPtr = nullptr
                ) const;

                //- Check face warpage: decompose face and check ratio between
                //  magnitude of sum of triangle areas and sum of magnitude of
                //  triangle areas.
                virtual bool checkFaceFlatness
                (
                    const bool report,
                    const scalar warnFlatness,  // When to include in set.
                    labelHashSet* setPtr
                ) const;

                //- Check for point-point-nearness,
                //  e.g. colocated points which may be part of baffles.
                virtual bool checkPointNearness
                (
                    const bool report,
                    const scalar reportDistSqr,
                    labelHashSet* setPtr = nullptr
                ) const;

                //- Check edge length
                virtual bool checkEdgeLength
                (
                    const bool report,
                    const scalar minLenSqr,
                    labelHashSet* setPtr = nullptr
                ) const;

                //- Check for concave cells by the planes of faces
                virtual bool checkConcaveCells
                (
                    const bool report = false,
                    labelHashSet* setPtr = nullptr
                ) const;


            //- Check mesh topology for correctness.
            //  Returns false for no error.
            virtual bool checkTopology(const bool report = false) const;

            //- Check mesh geometry (& implicitly topology) for correctness.
            //  Returns false for no error.
            virtual bool checkGeometry(const bool report = false) const;

            //- Check mesh for correctness. Returns false for no error.
            virtual bool checkMesh(const bool report = false) const;

            //- Set the closedness ratio warning threshold
            static scalar setClosedThreshold(const scalar);

            //- Set the aspect ratio warning threshold
            static scalar setAspectThreshold(const scalar);

            //- Set the non-orthogonality warning threshold in degrees
            static scalar setNonOrthThreshold(const scalar);

            //- Set the skewness warning threshold as percentage
            //  of the face area vector
            static scalar setSkewThreshold(const scalar);


        // Useful derived info

            //- Return true if the point in the cell bounding box.
            //  The bounding box may be isotropically inflated by the fraction
            //  inflationFraction
            bool pointInCellBB
            (
                const point& p,
                label celli,
                scalar inflationFraction = 0
            ) const;

            //- Return true if the  point is in the cell
            bool pointInCell(const point& p, label celli) const;

            //- Find the cell with the nearest cell centre to location
            label findNearestCell(const point& location) const;

            //- Find cell enclosing this location (-1 if not in mesh)
            label findCell(const point& location) const;


        //  Storage management

            //- Print a list of all the currently allocated mesh data
            void printAllocated() const;

            // Per storage whether allocated
            inline bool hasCellShapes() const;
            inline bool hasEdges() const;
            inline bool hasCellCells() const;
            inline bool hasEdgeCells() const;
            inline bool hasPointCells() const;
            inline bool hasCells() const;
            inline bool hasEdgeFaces() const;
            inline bool hasPointFaces() const;
            inline bool hasCellEdges() const;
            inline bool hasFaceEdges() const;
            inline bool hasPointEdges() const;
            inline bool hasPointPoints() const;
            inline bool hasCellPoints() const;
            inline bool hasCellCentres() const;
            inline bool hasFaceCentres() const;
            inline bool hasCellVolumes() const;
            inline bool hasFaceAreas() const;

            // On-the-fly addressing calculation. These functions return either
            // a reference to the full addressing (if already calculated) or
            // a reference to the supplied storage. The one-argument ones
            // use member DynamicList labels_ so be careful when not storing
            // result.

            //- cellCells using cells.
            const labelList& cellCells
            (
                const label celli,
                DynamicList<label>&
            ) const;

            const labelList& cellCells(const label celli) const;

            //- cellPoints using cells
            const labelList& cellPoints
            (
                const label celli,
                DynamicList<label>&
            ) const;

            const labelList& cellPoints(const label celli) const;

            //- pointCells using pointFaces
            const labelList& pointCells
            (
                const label pointi,
                DynamicList<label>&
            ) const;

            const labelList& pointCells(const label pointi) const;

            //- pointPoints using edges, pointEdges
            const labelList& pointPoints
            (
                const label pointi,
                DynamicList<label>&
            ) const;

            const labelList& pointPoints(const label pointi) const;

            //- faceEdges using pointFaces, edges, pointEdges
            const labelList& faceEdges
            (
                const label facei,
                DynamicList<label>&
            ) const;

            const labelList& faceEdges(const label facei) const;

            //- edgeFaces using pointFaces, edges, pointEdges
            const labelList& edgeFaces
            (
                const label edgeI,
                DynamicList<label>&
            ) const;

            const labelList& edgeFaces(const label edgeI) const;

            //- edgeCells using pointFaces, edges, pointEdges
            const labelList& edgeCells
            (
                const label edgeI,
                DynamicList<label>&
            ) const;

            const labelList& edgeCells(const label edgeI) const;

            //- cellEdges using cells, pointFaces, edges, pointEdges
            const labelList& cellEdges
            (
                const label celli,
                DynamicList<label>&
            ) const;

            const labelList& cellEdges(const label celli) const;


            //- Clear geometry
            void clearGeom();

            //- Clear topological data
            void clearAddressing();

            //- Clear all geometry and addressing unnecessary for CFD
            void clearOut();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const primitiveMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "primitiveMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
