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
    Foam::edgeCollapser

Description
    Does polyTopoChanges to remove edges. Can remove faces due to edge
    collapse but can not remove cells due to face removal!
    Also removes unused points.

SourceFiles
    edgeCollapser.C

\*---------------------------------------------------------------------------*/

#ifndef edgeCollapser_H
#define edgeCollapser_H

#include "pointEdgeCollapse.H"
#include "DynamicList.H"
#include "Field.H"
#include "pointFieldFwd.H"
#include "Map.H"
#include "labelPair.H"
#include "HashSet.H"
#include "typeInfo.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class PackedBoolList;
class polyTopoChange;
class globalIndex;
class face;
class edge;

/*---------------------------------------------------------------------------*\
                        Class edgeCollapser Declaration
\*---------------------------------------------------------------------------*/

class edgeCollapser
{
public:

    // The type of collapse of a face
    enum collapseType
    {
        noCollapse = 0,
        toPoint    = 1,
        toEdge     = 2
    };


private:

    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Controls collapse of a face to an edge
        const scalar guardFraction_;

        //- Only collapse face to a point if high aspect ratio
        const scalar maxCollapseFaceToPointSideLengthCoeff_;

        //- Allow a face to be collapsed to a point early, before the test
        //  to collapse to an edge
        const Switch allowEarlyCollapseToPoint_;

        //- Fraction of maxCollapseFaceToPointSideLengthCoeff_ to use when
        //  allowEarlyCollapseToPoint_ is on
        const scalar allowEarlyCollapseCoeff_;


    // Private Member Functions

        //- Create an edgeList of edges in facei which have both their points
        //  in pointLabels
        labelList edgesFromPoints
        (
            const label& facei,
            const labelList& pointLabels
        ) const;

        //- Collapse a face to an edge, marking the collapsed edges and new
        //  locations for points that will move as a result of the collapse
        void collapseToEdge
        (
            const label facei,
            const pointField& pts,
            const labelList& pointPriority,
            const vector& collapseAxis,
            const point& fC,
            const labelList& facePtsNeg,
            const labelList& facePtsPos,
            const scalarList& dNeg,
            const scalarList& dPos,
            const scalar dShift,
            PackedBoolList& collapseEdge,
            Map<point>& collapsePointToLocation
        ) const;

        //- Collapse a face to a point, marking the collapsed edges and new
        //  locations for points that will move as a result of the collapse
        void collapseToPoint
        (
            const label& facei,
            const pointField& pts,
            const labelList& pointPriority,
            const point& fC,
            const labelList& facePts,
            PackedBoolList& collapseEdge,
            Map<point>& collapsePointToLocation
        ) const;

        //- Do an eigenvector analysis of the face to get its collapse axis
        //  and aspect ratio
        void faceCollapseAxisAndAspectRatio
        (
            const face& f,
            const point& fC,
            vector& collapseAxis,
            scalar& aspectRatio
        ) const;

        //- Return the target length scale for each face
        scalarField calcTargetFaceSizes() const;

        //- Decides whether a face should be collapsed (and if so it it is to a
        //  point or an edge)
        collapseType collapseFace
        (
            const labelList& pointPriority,
            const face& f,
            const label facei,
            const scalar targetFaceSize,
            PackedBoolList& collapseEdge,
            Map<point>& collapsePointToLocation,
            const scalarField& faceFilterFactor
        ) const;

        //- Return label of point that has the highest priority. This will be
        //  the point on the edge that will be collapsed to.
        label edgeMaster(const labelList& pointPriority, const edge& e) const;

        //- Decides which points in an edge to collapse, based on their priority
        void checkBoundaryPointMergeEdges
        (
            const label pointi,
            const label otherPointi,
            const labelList& pointPriority,
            Map<point>& collapsePointToLocation
        ) const;

        //- Helper function that breaks strings of collapses if an edge is not
        //  labelled to collapse, but its points both collapse to the same
        //  location
        label breakStringsAtEdges
        (
            PackedBoolList& collapseEdge,
            List<pointEdgeCollapse>& allPointInfo
        ) const;

        //- Prevent face pinching by finding points in a face that will be
        //  collapsed to the same location, but that are not ordered
        //  consecutively in the face
        void determineDuplicatePointsOnFace
        (
            const face& f,
            PackedBoolList& markedPoints,
            labelHashSet& uniqueCollapses,
            labelHashSet& duplicateCollapses,
            List<pointEdgeCollapse>& allPointInfo
        ) const;

        //- Count the number of edges on the face that will exist as a result
        //  of the collapse
        label countEdgesOnFace
        (
            const face& f,
            List<pointEdgeCollapse>& allPointInfo
        ) const;

        //- Does the face have fewer than 3 edges as a result of the potential
        //  collapse
        bool isFaceCollapsed
        (
            const face& f,
            List<pointEdgeCollapse>& allPointInfo
        ) const;

        //- Given the collapse information, propagates the information using
        //  PointEdgeWave. Result is a list of new point locations and indices
        label syncCollapse
        (
            const globalIndex& globalPoints,
            const labelList& boundaryPoint,
            const PackedBoolList& collapseEdge,
            const Map<point>& collapsePointToLocation,
            List<pointEdgeCollapse>& allPointInfo
        ) const;

        //- Renumber f with new vertices. Removes consecutive duplicates
        void filterFace
        (
            const Map<DynamicList<label>>& collapseStrings,
            const List<pointEdgeCollapse>& allPointInfo,
            face& f
        ) const;


public:

    //- Runtime type information
    ClassName("edgeCollapser");


    // Constructors

        //- Construct from mesh
        edgeCollapser(const polyMesh& mesh);

        //- Construct from mesh and dict
        edgeCollapser(const polyMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        edgeCollapser(const edgeCollapser&) = delete;


    // Member Functions

        // Check

            //- Calls motionSmoother::checkMesh and returns a set of bad faces
            static HashSet<label> checkBadFaces
            (
                const polyMesh& mesh,
                const dictionary& meshQualityDict
            );

            //- Check mesh and mark points on faces in error
            //  Returns boolList with points in error set
            static label checkMeshQuality
            (
                const polyMesh& mesh,
                const dictionary& meshQualityDict,
                PackedBoolList& isErrorPoint
            );

            //- Ensure that the collapse is parallel consistent and update
            //  allPointInfo.
            //  Returns a list of edge collapses that is consistent across
            //  coupled boundaries and a list of pointEdgeCollapses.
            void consistentCollapse
            (
                const globalIndex& globalPoints,
                const labelList& pointPriority,
                const Map<point>& collapsePointToLocation,
                PackedBoolList& collapseEdge,
                List<pointEdgeCollapse>& allPointInfo,
                const bool allowCellCollapse = false
            ) const;


        // Query

            //- Play commands into polyTopoChange to create mesh.
            //  Return true if anything changed.
            bool setRefinement
            (
                const List<pointEdgeCollapse>& allPointInfo,
                polyTopoChange& meshMod
            ) const;

            //- Mark (in collapseEdge) any edges to collapse
            label markSmallEdges
            (
                const scalarField& minEdgeLen,
                const labelList& pointPriority,
                PackedBoolList& collapseEdge,
                Map<point>& collapsePointToLocation
            ) const;

            //- Mark (in collapseEdge) any edges to merge
            label markMergeEdges
            (
                const scalar maxCos,
                const labelList& pointPriority,
                PackedBoolList& collapseEdge,
                Map<point>& collapsePointToLocation
            ) const;

            //- Find small faces and sliver faces in the mesh and mark the
            //  edges that need to be collapsed in order to remove these faces.
            //  Also returns a map of new locations for points that will move
            //  as a result of the collapse.
            //  Use in conjunctions with edgeCollapser to synchronise the
            //  collapses and modify the mesh
            labelPair markSmallSliverFaces
            (
                const scalarField& faceFilterFactor,
                const labelList& pointPriority,
                PackedBoolList& collapseEdge,
                Map<point>& collapsePointToLocation
            ) const;

            //- Marks edges in the faceZone indirectPatchFaces for collapse
            labelPair markFaceZoneEdges
            (
                const faceZone& fZone,
                const scalarField& faceFilterFactor,
                const labelList& pointPriority,
                PackedBoolList& collapseEdge,
                Map<point>& collapsePointToLocation
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const edgeCollapser&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
