/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::TriPatchIntersection

Description
    Patch intersection based on triangular faces. Intersects and combines two
    triangulated patches incrementally. The intersected surface is valid at
    every stage of the process. Failure to intersect does not produce a
    catastrophic error. Rather, it results in regions of the surface remaining
    associated with only one of the source or the target patch.

SourceFiles
    TriPatchIntersection.C

\*---------------------------------------------------------------------------*/

#ifndef TriPatchIntersection_H
#define TriPatchIntersection_H

#include "PatchIntersection.H"
#include "star.H"
#include "polygonTriangulate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class TriPatchIntersection Declaration
\*---------------------------------------------------------------------------*/

template<class SrcPatchType, class TgtPatchType>
class TriPatchIntersection
:
    public PatchIntersection<SrcPatchType, TgtPatchType>
{
private:

    // Private Data

        // Points

            //- The source points
            DynamicField<point> srcPoints_;

            //- The source point normals
            DynamicField<vector> srcPointNormals_;

            //- The target points. Reference to base class data.
            DynamicField<point>& tgtPoints_;

            //- Point-point addressing. Facilitates quick point removal.
            //  Instead of shuffling up, when points are combined this just
            //  acts as a redirect. Everything gets resolved on clean.
            DynamicList<label> pointPoints_;


        // Edges

            //- The two triangles adjacent to each edge
            DynamicList<labelPair> edgeTris_;

            //- ...
            DynamicList<labelPair> intersectEdgeFaces_;

            //- ...
            DynamicList<labelPair> nonIntersectEdgeFaces_;


        // Tris

            //- The triangles' points
            DynamicList<triFace> triPoints_;

            //- The triangles' edges
            DynamicList<FixedList<label, 3>> triEdges_;

            //- The source patch face associated with each triangle. -1 if not
            //  associated with a source face.
            DynamicList<label> triSrcFace_;

            //- The target patch face associated with each triangle. -1 if not
            //  associated with a target face.
            DynamicList<label> triTgtFace_;

            //- The triangles constructed from each source patch face
            List<DynamicList<label>> srcFaceTris_;

            //- The triangles constructed from each target patch face
            List<DynamicList<label>> tgtFaceTris_;


        // Faces

            //- ...
            DynamicList<labelList> faceEdges_;


        // Removal

            //- Indices of triangles that have been removed. Prevents shuffling
            //  up and re-allocation during intersection. Cleared out on clean.
            DynamicList<label> removedTris_;

            //- Indices of edges that have been removed. Prevents shuffling
            //  up and re-allocation during intersection. Cleared out on clean.
            DynamicList<label> removedEdges_;


        // Front propagation

            //- The edges in the front along which the intersection propagates
            DynamicList<label> frontEdgeEdges_;

            //- Inverse of the above
            DynamicList<label> edgeFrontEdges_;


        // Insertion

            //- Indices of triangles considered candidates for insertion
            DynamicList<label> candidateTriTris_;

            //- Inverse of the above
            DynamicList<label> triCandidateTris_;


        // Marking

            //- Indices of triangles that have been "marked". Used for multiple
            //  purposes.
            DynamicList<label> markedTriTris_;

            //- Inverse of the above
            DynamicList<label> triMarkedTris_;


        // Triangulation

            //- Triangulation engine
            polygonTriangulate polygonTriangulate_;


        // Polygon generation

            //- Star propagation engine
            star star_;


        // Debugging

            //- Iteration counter with which to number surface files
            label writei_;


    // Private Member Functions

        //- Check consistency of the mesh
        void checkPatchFace(const label patchFacei, const bool isSrc) const;
        void checkPatchEdge(const label patchFacei, const bool isSrc) const;
        void checkPatchFaces(const bool isSrc) const;

        //- Remove an edge
        void removeEdge(const label edgei);

        //- Remove a tri
        void removeTri(const label trii);

        //- Create a new triangle and return it's index
        label newTrii();

        //- Create a new edge and return it's index
        label newEdgei();

        //- Create a new point and return it's index
        label newPointi();

        //- Return a point of a given triangle
        label triPoint(const label trii, const label triPointi) const;

        //- Return the points of a given triangle
        triFace triPoints(const label trii) const;

        //- Get the values on the points of a given triangle
        template <class Type>
        FixedList<Type, 3> triPointValues
        (
            const label trii,
            const UList<Type> values
        ) const;

        //- Get a list indicating whether a tri "owns" it's edges; i.e.,
        //  whether they are numbered in the same order as the points
        FixedList<bool, 3> triOwns(const label trii) const;

        //- Get the shared point indices of another tri in a tri, or -1 if
        //  the point is not shared
        FixedList<label, 3> triOtherTriPoints
        (
            const label trii,
            const label otherTrii
        ) const;

        //- Get the shared edge indices of another tri in a tri, or -1 if
        //  the point is not shared
        FixedList<label, 3> triOtherTriEdges
        (
            const label trii,
            const label otherTrii
        ) const;

        //- Return the tri-edge-points for a given tri-edge
        edge triEdgePoints(const label trii, const label triEdgei) const;

        //- Return the edge-points for a given edge
        edge edgePoints(const label edgei) const;

        //- Return a point for the given patch face
        label patchFacePoint
        (
            const label patchFacei,
            const label patchFacePointi,
            const bool isSrc
        ) const;

        //- Return the points for the given patch face
        triFace patchFacePoints
        (
            const label patchFacei,
            const bool isSrc
        ) const;

        //- Get the values on the points of a given patch face
        template <class Type>
        FixedList<Type, 3> patchFacePointValues
        (
            const label patchFacei,
            const bool isSrc,
            const UList<Type>& values
        ) const;

        //- Get a list indicating whether a patch face "owns" it's edges; i.e.,
        //  whether they are numbered in the same order as the points
        FixedList<bool, 3> patchFaceOwns
        (
            const label patchFacei,
            const bool isSrc
        ) const;

        //- Get the shared point indices of another patch face in a patch face,
        //  or -1 if the point is not shared
        FixedList<label, 3> patchFaceOtherPatchFacePoints
        (
            const label patchFacei,
            const label otherPatchFacei,
            const bool isSrc
        ) const;

        //- Return a patch point for the given patch face
        label patchFacePatchPoint
        (
            const label patchFacei,
            const label patchFacePointi,
            const bool isSrc
        ) const;

        //- Return the patch points for the given patch face
        triFace patchFacePatchPoints
        (
            const label patchFacei,
            const bool isSrc
        ) const;

        //- Return a patch edge for the given patch face
        label patchFacePatchEdge
        (
            const label patchFacei,
            const label patchFaceEdgei,
            const bool isSrc
        ) const;

        //- Return the patch edges for the given patch face
        triFace patchFacePatchEdges
        (
            const label patchFacei,
            const bool isSrc
        ) const;

        //- Compute the patch edge that a given edge lies along
        label edgePatchEdge(const label edgei, const bool isSrc) const;

        //- Compute the patch edges that a given edge lies along
        labelPair edgePatchEdges(const label edgei) const;

        //- Add a tri. Return the new tri label.
        label addTri
        (
            const triFace& pointis,
            const FixedList<label, 3>& edgeis,
            const label patchFacei,
            const bool isSrc
        );

        //- Flip an edge
        void flipEdge(const label edgei);

        //- Return the signed distance squared of the given point outside of a
        //  triangle's circumcircle. A negative value means the point is inside
        //  the circle.
        scalar circumDistSqr(const label trii, const label pointi) const;

        //- Insert points into the given tri or edge. If an edge, points are to
        //  be given in order along the edge.
        void insertPoints
        (
            const label triOrEdgei,
            const bool isTri,
            const UList<label>& pointis,
            UList<label>& insertionEdgeis,
            const UList<label>& fixedEdgeis
        );

        //- Return whether or not a point can be intersected with the other side
        bool pointCanIntersect(const label pointi) const;

        //- Return whether or not an edge can be intersected with the other side
        bool edgeCanIntersect(const label pointi) const;

        /*
        //- Snap points to points and edges between two tris
        //  !!! Not a good idea. Gets tangled up. Need to snap entire faces.
        bool snapTris(const label srcTrii, const label tgtTrii);
        */

        //- Snap points to points and edges between two patch faces
        void snapPatchFaceTris
        (
            const label srcFacei,
            const label tgtFacei,
            const scalar snapTol
        );

        //- Intersect two tris
        bool intersectTris(const label srcTrii, const label tgtTrii);

        //- Intersect the triangulations of the given patch faces
        void intersectPatchFaceTris
        (
            const label srcFacei,
            const label tgtFacei
        );

        //- Make the triangulations of the given patch faces conform
        bool conformPatchFaceTris
        (
            const label patchFacei,
            const label otherPatchFacei,
            const bool isSrc
        );

        //- Make the triangulations of the given patch faces conform
        bool conformPatchFaceTris
        (
            const label srcFacei,
            const label tgtFacei
        );

        //- Combine conformal parts of the given patch faces into intersection
        //  faces
        bool combinePatchFaceTris(const label srcFacei, const label tgtFacei);

        //- Initialise the member data
        void initialise(const vectorField& srcPointNormals);

        //- Clean up removed tris and edges
        void clean();

        //- Finalise the data in the base class
        void finalise();

        //- Undo the finalise process so the intersection process can continue
        void unFinalise();

        //- Write the current state of the surfaces for debugging purposes
        void write();

        //- Write the current state of the given patch face triangulation for
        //  debugging purposes
        void writePatchFace(const label patchFacei, const bool isSrc) const;


public:

    // Runtime type information

        virtual word type() const
        {
            return "tri" + patchIntersection::typeName.capitalise();
        }


    // Constructors

        //- Construct from a source and a target patch
        TriPatchIntersection
        (
            const SrcPatchType& srcPatch,
            const TgtPatchType& tgtPatch,
            const scalar snapTol
        );

        //- Construct from a source and a target patch, and specified source
        //  point normals
        TriPatchIntersection
        (
            const SrcPatchType& srcPatch,
            const vectorField& srcPointNormals,
            const TgtPatchType& tgtPatch,
            const scalar snapTol
        );


    // Member Functions

        // Access

            //- ...
            inline const DynamicList<labelList>& faceEdges() const
            {
                return faceEdges_;
            }

            //- ...
            inline const DynamicList<labelPair>& intersectEdgeFaces() const
            {
                return intersectEdgeFaces_;
            }

            //- ...
            inline const DynamicList<labelPair>& nonIntersectEdgeFaces() const
            {
                return nonIntersectEdgeFaces_;
            }


    //- Destructor
    virtual ~TriPatchIntersection();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TriPatchIntersection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
