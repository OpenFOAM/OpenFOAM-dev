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
    Foam::PatchTools

Description
    A collection of tools for searching, sorting PrimitivePatch information.

    The class could also be extended to include more that just static methods.

SourceFiles
    PatchTools.C
    PatchToolsCheck.C
    PatchToolsEdgeOwner.C
    PatchToolsGatherAndMerge.C
    PatchToolsMatch.C
    PatchToolsNormals.C
    PatchToolsSearch.C
    PatchToolsSortEdges.C
    PatchToolsSortPoints.C

\*---------------------------------------------------------------------------*/

#ifndef PatchTools_H
#define PatchTools_H

#include "PrimitivePatch.H"
#include "HashSet.H"
#include "globalIndex.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class PackedBoolList;
class boundBox;

/*---------------------------------------------------------------------------*\
                         Class PatchTools Declaration
\*---------------------------------------------------------------------------*/

class PatchTools
{
public:

    //- Check for orientation issues.
    //  Returns true if problems were found.
    //  If a normal flips across an edge, places it in the HashSet
    template<class FaceList, class PointField>
    static bool checkOrientation
    (
        const PrimitivePatch<FaceList, PointField>&,
        const bool report = false,
        labelHashSet* marked = 0
    );


    //- Fill faceZone with currentZone for every face reachable
    //  from facei without crossing edge marked in borderEdge.
    //  Note: faceZone has to be sized nFaces before calling.
    template<class BoolListType, class FaceList, class PointField>
    static void markZone
    (
        const PrimitivePatch<FaceList, PointField>&,
        const BoolListType& borderEdge,
        const label facei,
        const label currentZone,
        labelList& faceZone
    );

    //- Size and fills faceZone with zone of face.
    //  Zone is area reachable by edge crossing without crossing borderEdge.
    //  Returns number of zones.
    template<class BoolListType, class FaceList, class PointField>
    static label markZones
    (
        const PrimitivePatch<FaceList, PointField>&,
        const BoolListType& borderEdge,
        labelList& faceZone
    );

    //- Determine the mapping for a sub-patch.
    //  Only include faces for which bool-list entry is true.
    //  \param[in]  p patch to be searched on
    //  \param[in]  includeFaces faces to include
    //  \param[out] pointMap mapping new to old localPoints
    //  \param[out] faceMap  mapping new to old faces
    template<class BoolListType, class FaceList, class PointField>
    static void subsetMap
    (
        const PrimitivePatch<FaceList, PointField>& p,
        const BoolListType& includeFaces,
        labelList& pointMap,
        labelList& faceMap
    );

    //-
    template<class FaceList, class PointField>
    static void calcBounds
    (
        const PrimitivePatch<FaceList, PointField>& p,
        boundBox& bb,
        label& nPoints
    );

    //- Return edge-face addressing sorted by angle around the edge.
    //  Orientation is anticlockwise looking from edge.vec(localPoints())
    template<class FaceList, class PointField>
    static labelListList sortedEdgeFaces
    (
        const PrimitivePatch<FaceList, PointField>&
    );

    //- Return point-edge addressing sorted by order around the point.
    template<class FaceList, class PointField>
    static labelListList sortedPointEdges
    (
        const PrimitivePatch<FaceList, PointField>&
    );

    //- If 2 face neighbours: label of face where ordering of edge
    //  is consistent with righthand walk.
    //  If 1 neighbour: label of only face.
    //  If >2 neighbours: undetermined.
    template<class FaceList, class PointField>
    static labelList edgeOwner
    (
        const PrimitivePatch<FaceList, PointField>&
    );


    //- Find corresponding points on patches sharing the same points
    //  p1PointLabels : points on p1 that were matched
    //  p2PointLabels : corresponding points on p2
    template
    <
        class FaceList1,
        class PointField1,
        class FaceList2,
        class PointField2
    >
    static void matchPoints
    (
        const PrimitivePatch<FaceList1, PointField1>& p1,
        const PrimitivePatch<FaceList2, PointField2>& p2,

        labelList& p1PointLabels,
        labelList& p2PointLabels
    );

    //- Find corresponding edges on patches sharing the same points
    //  p1EdgeLabels    : edges on p1 that were matched
    //  p2EdgeLabels    : corresponding edges on p2
    //  sameOrientation : same orientation?
    template
    <
        class FaceList1,
        class PointField1,
        class FaceList2,
        class PointField2
    >
    static void matchEdges
    (
        const PrimitivePatch<FaceList1, PointField1>& p1,
        const PrimitivePatch<FaceList2, PointField2>& p2,

        labelList& p1EdgeLabels,
        labelList& p2EdgeLabels,
        PackedBoolList& sameOrientation
    );


    //- Return parallel consistent point normals for patches using mesh points.
    template<class FaceList, class PointField>
    static tmp<pointField> pointNormals
    (
        const polyMesh&,
        const PrimitivePatch<FaceList, PointField>&
    );


    //- Return parallel consistent edge normals for patches using mesh points.
    //  Supply with patch matching info from matchEdges.
    template<class FaceList, class PointField>
    static tmp<pointField> edgeNormals
    (
        const polyMesh&,
        const PrimitivePatch<FaceList, PointField>&,
        const labelList& patchEdges,
        const labelList& coupledEdges
    );


    //- Gather points and faces onto master and merge into single patch.
    //  Note: uses faces/points, not localFaces/localPoints.
    template<class FaceList, class PointField>
    static void gatherAndMerge
    (
        const scalar mergeDist,
        const PrimitivePatch<FaceList, PointField>& p,
        Field<typename PrimitivePatch<FaceList, PointField>::PointType>&
            mergedPoints,
        List<typename PrimitivePatch<FaceList, PointField>::FaceType>&
            mergedFaces,
        labelList& pointMergeMap
    );

    //- Gather (mesh!) points and faces onto master and merge collocated
    //  points into a single patch. Uses coupled point mesh
    //  structure so does not need tolerances.
    //  On master and slave returns:
    //  - pointToGlobal : for every local point index the global point index
    //  - uniqueMeshPointLabels : my local mesh points
    //  - globalPoints : global numbering for the global points
    //  - globalFaces : global numbering for the faces
    //  On master only:
    //  - mergedFaces : the merged faces
    //  - mergedPoints : the merged points
    template<class FaceList>
    static void gatherAndMerge
    (
        const polyMesh& mesh,
        const FaceList& faces,
        const labelList& meshPoints,
        const Map<label>& meshPointMap,

        labelList& pointToGlobal,
        labelList& uniqueMeshPointLabels,
        autoPtr<globalIndex>& globalPoints,
        autoPtr<globalIndex>& globalFaces,
        List<typename FaceList::value_type>& mergedFaces,
        pointField& mergedPoints
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PatchTools.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
