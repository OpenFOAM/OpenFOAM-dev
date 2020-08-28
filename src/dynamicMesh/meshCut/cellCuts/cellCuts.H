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
    Foam::cellCuts

Description
    Description of cuts across cells.

    Description of cut is given as list of vertices and
    list of edges to be cut (and position on edge).

    Does some checking of correctness/non-overlapping of cuts. 2x2x2
    refinement has to be done in three passes since cuts can not overlap
    (would make addressing too complicated)

    Introduces concept of 'cut' which is either an existing vertex
    or a edge.

    Input can either be
    -# list of cut vertices and list of cut edges. Constructs cell
       circumference walks ('cellLoops').

    -# list of cell circumference walks. Will filter them so they
       don't overlap.

    -# cellWalker and list of cells to refine (refineCell). Constructs
       cellLoops and does B. cellWalker is class which can cut a single
       cell using a plane through the cell centre and in a certain normal
       direction

    CellCuts constructed from cellLoops (B, C) can have multiple cut-edges
    and/or cut-point as long as there is per face only one (or none) cut across
    a face, i.e. a face can only be split into two faces.

    The information available after construction:
      - pointIsCut, edgeIsCut.
      - faceSplitCut : the cross-cut of a face.
      - cellLoops : the loop which cuts across a cell
      - cellAnchorPoints : per cell the vertices on one side which are
        considered the anchor points.

    AnchorPoints: connected loops have to be oriented in the same way to
    be able to grow new internal faces out of the same bottom faces.
    (limitation of the mapping procedure). The loop is cellLoop is oriented
    such that the normal of it points towards the anchorPoints.
    This is the only routine which uses geometric information.


    Limitations:
    - cut description is very precise. Hard to get right.
    - do anchor points need to be unique per cell? Very inefficient.
    - is orientation guaranteed to be correct? Uses geometric info so can go
      wrong on highly distorted meshes.
    - is memory inefficient. Probably need to use Maps instead of
      labelLists.

SourceFiles
    cellCuts.C

\*---------------------------------------------------------------------------*/

#ifndef cellCuts_H
#define cellCuts_H

#include "edgeVertex.H"
#include "labelList.H"
#include "boolList.H"
#include "scalarField.H"
#include "pointField.H"
#include "DynamicList.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class cellLooper;
class refineCell;
class plane;

/*---------------------------------------------------------------------------*\
                          Class cellCuts Declaration
\*---------------------------------------------------------------------------*/

class cellCuts
:
    public edgeVertex
{
    // Private Data

        // Per point/edge status

            //- Is mesh point cut
            boolList pointIsCut_;

            //- Is edge cut
            boolList edgeIsCut_;

            //- If edge is cut gives weight (0->start() to 1->end())
            scalarField edgeWeight_;


        // Cut addressing

            //- Cuts per existing face (includes those along edge of face)
            //  Cuts in no particular order.
            mutable autoPtr<labelListList> faceCutsPtr_;

            //- Per face : cut across edge (so not along existing edge)
            //  (can only be one per face)
            Map<edge> faceSplitCut_;


        // Cell-cut addressing

            //- Loop across cell circumference
            labelListList cellLoops_;

            //- Number of valid loops in cellLoops_
            label nLoops_;

            //- For each cut cell the points on the 'anchor' side of the cell.
            labelListList cellAnchorPoints_;


    // Private Static Functions

        //- Find value in first n elements of list.
        static label findPartIndex
        (
            const labelList&,
            const label n,
            const label val
        );

        //- Create boolList with all labels specified set to true
        //  (and rest to false)
        static boolList expand(const label size, const labelList& labels);

        //- Create scalarField with all specified labels set to corresponding
        //  value in scalarField.
        static scalarField expand
        (
            const label,
            const labelList&,
            const scalarField&
        );

        //- Returns -1 or index of first element of lst that cannot be found
        //  in map.
        static label firstUnique
        (
            const labelList& lst,
            const Map<label>&
        );

    // Private Member Functions

        //- Debugging: write cell's edges and any cut vertices and edges
        //  (so no cell loop determined yet)
        void writeUncutOBJ(const fileName&, const label celli) const;

        //- Debugging: write cell's edges, loop and anchors to directory.
        void writeOBJ
        (
            const fileName& dir,
            const label celli,
            const pointField& loopPoints,
            const labelList& anchors
        ) const;

        //- Find face on cell using the two edges.
        label edgeEdgeToFace
        (
            const label celli,
            const label edgeA,
            const label edgeB
        ) const;


        //- Find face on cell using an edge and a vertex.
        label edgeVertexToFace
        (
            const label celli,
            const label edgeI,
            const label vertI
        ) const;

        //- Find face using two vertices (guaranteed not to be along edge)
        label vertexVertexToFace
        (
            const label celli,
            const label vertA,
            const label vertB
        ) const;


        // Cut addressing

            //- Calculate faceCuts in face vertex order.
            void calcFaceCuts() const;


        // Loop (cuts on cell circumference) calculation

            //- Find edge (or -1) on facei using vertices v0,v1
            label findEdge
            (
                const label facei,
                const label v0,
                const label v1
            ) const;

            //- Find face on which all cuts are (very rare) or -1.
            label loopFace(const label celli, const labelList& loop) const;

            //- Cross otherCut into next faces (not exclude0, exclude1)
            bool walkPoint
            (
                const label celli,
                const label startCut,

                const label exclude0,
                const label exclude1,

                const label otherCut,

                label& nVisited,
                labelList& visited
            ) const;

            //- Cross cut (which is edge on facei) onto next face
            bool crossEdge
            (
                const label celli,
                const label startCut,
                const label facei,
                const label otherCut,

                label& nVisited,
                labelList& visited
            ) const;

            // wrapper around visited[nVisited++] = cut. Checks for duplicate
            // cuts.
            bool addCut
            (
                const label celli,
                const label cut,
                label& nVisited,
                labelList& visited
            ) const;

            //- Walk across facei following cuts, starting at cut. Stores cuts
            //  visited
            // Returns true if valid walk.
            bool walkFace
            (
                const label celli,
                const label startCut,
                const label facei,
                const label cut,

                label& lastCut,
                label& beforeLastCut,
                label& nVisited,
                labelList& visited
            ) const;

            //- Walk across cuts (cut edges or cut vertices) of cell. Stops when
            //  hit cut  already visited. Returns true when loop of 3 or more
            //  vertices found.
            bool walkCell
            (
                const label celli,
                const label startCut,   // overall starting cut
                const label facei,
                const label prevCut,    // current cut
                label& nVisited,
                labelList& visited
            ) const;

            //- Determine for every cut cell the face it is cut by.
            void calcCellLoops(const labelList& cutCells);


        // Cell anchoring

            //- Are there enough faces on anchor side of celli?
            bool checkFaces
            (
                const label celli,
                const labelList& anchorPoints
            ) const;

            //- Walk unset edges of single cell from starting point and
            //  marks visited edges and vertices with status.
            void walkEdges
            (
                const label celli,
                const label pointi,
                const label status,

                Map<label>& edgeStatus,
                Map<label>& pointStatus
            ) const;

            //- Check anchor points on 'outside' of loop
            bool loopAnchorConsistent
            (
                const label celli,
                const pointField& loopPts,
                const labelList& anchorPoints
            ) const;

            //- Determines set of anchor points given a loop. The loop should
            //  split the cell into two. Returns true if a valid set of anchor
            //  points determined, false otherwise.
            bool calcAnchors
            (
                const label celli,
                const labelList& loop,
                const pointField& loopPts,

                labelList& anchorPoints
            ) const;

            //- Returns coordinates of points on loop with explicitly provided
            //  weights.
            pointField loopPoints
            (
                const labelList& loop,
                const scalarField& loopWeights
            ) const;

            //- Returns weights of loop. Inverse of loopPoints.
            scalarField loopWeights(const labelList& loop) const;

            //- Check if cut edges in loop are compatible with ones in
            //  edgeIsCut_
            bool validEdgeLoop
            (
                const labelList& loop,
                const scalarField& loopWeights
            ) const;

            //- Counts number of cuts on face.
            label countFaceCuts
            (
                const label facei,
                const labelList& loop
            ) const;

            //- Determine compatibility of loop with existing cut pattern.
            //  Does not use cut-addressing (faceCuts_, cutCuts_)
            bool conservativeValidLoop
            (
                const label celli,
                const labelList& loop
            ) const;

            //- Check if loop is compatible with existing cut pattern in
            //  pointIsCut, edgeIsCut, faceSplitCut.
            //  Calculates and returns for current cell the cut faces and the
            //  points on one side of the loop.
            bool validLoop
            (
                const label celli,
                const labelList& loop,
                const scalarField& loopWeights,
                Map<edge>& newFaceSplitCut,
                labelList& anchorPoints
            ) const;

            //- Update basic cut information from cellLoops.
            //  Assumes cellLoops_ and edgeWeight_ already set and consistent.
            //  Does not use any other information.
            void setFromCellLoops();

            //- Update basic cut information for single cell from cellLoop.
            bool setFromCellLoop
            (
                const label celli,
                const labelList& loop,
                const scalarField& loopWeights
            );

            //- Update basic cut information from cellLoops. Checks for
            //  consistency with existing cut pattern.
            void setFromCellLoops
            (
                const labelList& cellLabels,
                const labelListList& cellLoops,
                const List<scalarField>& cellLoopWeights
            );

            //- Update basic cut information to be consistent across
            //  coupled boundaries.
            void syncProc();

            //- Cut cells and update basic cut information from cellLoops.
            //  Checks each loop for consistency with existing cut pattern.
            void setFromCellCutter
            (
                const cellLooper&,
                const List<refineCell>& refCells
            );

            //- Same as above but now cut with prescribed plane (instead of
            //  just normal).
            void setFromCellCutter
            (
                const cellLooper&,
                const labelList& cellLabels,
                const List<plane>&
            );

            //- Set orientation of loops
            void orientPlanesAndLoops();

            //- Top level driver: addressing calculation and loop detection
            //  (loops splitting cells).
            void calcLoopsAndAddressing(const labelList& cutCells);

            //- Check various consistencies.
            void check() const;


public:

    //- Runtime type information
    ClassName("cellCuts");


    // Constructors

        //- Construct from cells to cut and pattern of cuts
        cellCuts
        (
            const polyMesh& mesh,
            const labelList& cutCells,
            const labelList& meshVerts,
            const labelList& meshEdges,
            const scalarField& meshEdgeWeights
        );

        //- Construct from pattern of cuts. Detect cells to cut.
        cellCuts
        (
            const polyMesh& mesh,
            const labelList& meshVerts,
            const labelList& meshEdges,
            const scalarField& meshEdgeWeights
        );

        //- Construct from complete cellLoops through specified cells.
        //  Checks for consistency.
        //  Constructs cut-cut addressing and cellAnchorPoints.
        cellCuts
        (
            const polyMesh& mesh,
            const labelList& cellLabels,
            const labelListList& cellLoops,
            const List<scalarField>& cellEdgeWeights
        );

        //- Construct from list of cells to cut and direction to cut in
        //  (always through cell centre) and celllooper.
        cellCuts
        (
            const polyMesh& mesh,
            const cellLooper& cellCutter,
            const List<refineCell>& refCells
        );

        //- Construct from list of cells to cut and plane to cut with and
        //  celllooper. (constructor above always cuts through cell centre)
        cellCuts
        (
            const polyMesh& mesh,
            const cellLooper& cellCutter,
            const labelList& cellLabels,
            const List<plane>& cutPlanes
        );

        //- Construct from components
        cellCuts
        (
            const polyMesh& mesh,
            const boolList& pointIsCut,
            const boolList& edgeIsCut,
            const scalarField& edgeWeight,
            const Map<edge>& faceSplitCut,
            const labelListList& cellLoops,
            const label nLoops,
            const labelListList& cellAnchorPoints
        );

        //- Disallow default bitwise copy construction
        cellCuts(const cellCuts&) = delete;


    //- Destructor
    ~cellCuts();

        //- Clear out demand driven storage
        void clearOut();


    // Member Functions

        // Access

            //- Is mesh point cut
            const boolList& pointIsCut() const
            {
                return pointIsCut_;
            }

            //- Is edge cut
            const boolList& edgeIsCut() const
            {
                return edgeIsCut_;
            }

            //- If edge is cut gives weight (ratio between start() and end())
            const scalarField& edgeWeight() const
            {
                return edgeWeight_;
            }

            //- Cuts per existing face (includes those along edge of face)
            //  Cuts in no particular order
            const labelListList& faceCuts() const
            {
                if (!faceCutsPtr_.valid())
                {
                    calcFaceCuts();
                }
                return faceCutsPtr_();
            }

            //- Gives for split face the two cuts that split the face into two.
            const Map<edge>& faceSplitCut() const
            {
                return faceSplitCut_;
            }

            //- For each cut cell the cut along the circumference.
            const labelListList& cellLoops() const
            {
                return cellLoops_;
            }

            //- Number of valid cell loops
            label nLoops() const
            {
                return nLoops_;
            }

            //- For each cut cell the points on the 'anchor' side of the cell.
            const labelListList& cellAnchorPoints() const
            {
                return cellAnchorPoints_;
            }


        // Other

            //- Returns coordinates of points on loop for given cell.
            //  Uses cellLoops_ and edgeWeight_
            pointField loopPoints(const label celli) const;

            //- Invert anchor point selection.
            labelList nonAnchorPoints
            (
                const labelList& cellPoints,
                const labelList& anchorPoints,
                const labelList& loop
            ) const;

            //- Flip loop for celli. Updates anchor points as well.
            void flip(const label celli);

            //- Flip loop for celli. Does not update anchors. Use with care
            //  (only if you're sure loop orientation is wrong)
            void flipLoopOnly(const label celli);


        // Write

            //- debugging:Write list of cuts to stream in OBJ format
            void writeOBJ
            (
                Ostream& os,
                const pointField& loopPoints,
                label& vertI
            ) const;

            //- debugging:Write all of cuts to stream in OBJ format
            void writeOBJ(Ostream& os) const;

            //- debugging:Write edges of cell and loop
            void writeCellOBJ(const fileName& dir, const label celli) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cellCuts&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
