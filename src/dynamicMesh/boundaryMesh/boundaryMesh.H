/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::boundaryMesh

Description
    Addressing for all faces on surface of mesh. Can either be read
    from polyMesh or from triSurface. Used for repatching existing meshes.

SourceFiles
    boundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryMesh_H
#define boundaryMesh_H

#include "bMesh.H"
#include "boundaryPatch.H"
#include "PrimitivePatch.H"
#include "PtrList.H"
#include "polyPatchList.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
class polyMesh;
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                           Class boundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class boundaryMesh
{
    // Static data

        //- Normal along which to divide faces into categories
        //  (used in getNearest)
        static const vector splitNormal_;

        //- Distance to face tolerance for getNearest. Triangles are considered
        //  near if they are nearer than distanceTol_*typDim where typDim is
        //  the largest distance from face centre to one of its vertices.
        static const scalar distanceTol_;

    // Private data

        //- All boundary mesh data. Reconstructed every time faces are repatched
        bMesh* meshPtr_;

        //- Patches. Reconstructed every time faces are repatched.
        PtrList<boundaryPatch> patches_;

        //- For every face in mesh() gives corresponding polyMesh face
        // (not sensible if mesh read from triSurface)
        labelList meshFace_;


        //
        // Edge handling
        //

        //- points referenced by feature edges.
        pointField featurePoints_;

        //- feature edges. Indices into featurePoints.
        edgeList featureEdges_;

        //- from feature edge to mesh edge.
        labelList featureToEdge_;

        //- from mesh edges to featureEdges_;
        labelList edgeToFeature_;

        //- Feature 'segments'. Collections of connected featureEdges.
        //  Indices into featureEdges_.
        labelListList featureSegments_;

        //- Additional edges (indices of mesh edges)
        labelList extraEdges_;


    // Private Member Functions

        //- Number of connected feature edges.
        label nFeatureEdges(label pointI) const;

        //- Step to next feature edge
        label nextFeatureEdge(const label edgeI, const label vertI) const;

        //- Return connected list of feature edges.
        labelList collectSegment
        (
            const boolList& isFeaturePoint,
            const label startEdgeI,
            boolList& featVisited
        ) const;

        //- Do point-edge walk to determine nearest (to edgeI). Stops if
        //  distance >= maxDistance. Used to determine edges close to seed
        //  point.
        void markEdges
        (
            const label maxDistance,
            const label edgeI,
            const label distance,
            labelList& minDistance,
            DynamicList<label>& visited
        ) const;

        //- Get index of polypatch by name
        label findPatchID(const polyPatchList&, const word&) const;

        //- Get index of patch for face
        label whichPatch(const polyPatchList&, const label) const;

        //- Gets labels of changed faces and propagates them to the edges.
        //  Returns labels of edges changed. Fills edgeRegion of visited edges
        //  with current region.
        labelList faceToEdge
        (
            const boolList& regionEdge,
            const label region,
            const labelList& changedFaces,
            labelList& edgeRegion
        ) const;

        //- Reverse of faceToEdge: gets edges and returns faces
        labelList edgeToFace
        (
            const label region,
            const labelList& changedEdges,
            labelList& faceRegion
        ) const;

        //- Finds area, starting at faceI, delimited by borderEdge. Marks all
        //  faces thus visited with currentZone.
        void markZone
        (
            const boolList& borderEdge,
            label faceI,
            label currentZone,
            labelList& faceZone
        ) const;


        //- Disallow default bitwise copy construct
        boundaryMesh(const boundaryMesh&);

        //- Disallow default bitwise assignment
        void operator=(const boundaryMesh&);


public:

    //- Runtime type information
    ClassName("boundaryMesh");


    // Constructors

        //- Construct null
        boundaryMesh();


    //- Destructor
    ~boundaryMesh();

        void clearOut();


    // Member Functions

        // Access

            const bMesh& mesh() const
            {
                if (!meshPtr_)
                {
                    FatalErrorIn("boundaryMesh::mesh()")
                        << "No mesh available. Probably mesh not yet"
                        << " read." << abort(FatalError);
                }
                return *meshPtr_;
            }

            const PtrList<boundaryPatch>& patches() const
            {
                return patches_;
            }


            //- Label of original face in polyMesh (before patchify(...))
            const labelList& meshFace() const
            {
                return meshFace_;
            }

            //- Feature points.
            const pointField& featurePoints() const
            {
                return featurePoints_;
            }

            //- Feature edges. Indices into featurePoints.
            const edgeList& featureEdges() const
            {
                return featureEdges_;
            }

            //- From index into featureEdge to index into meshedges,
            const labelList& featureToEdge() const
            {
                return featureToEdge_;
            }

            //- From edge into featureEdges
            const labelList& edgeToFeature() const
            {
                return edgeToFeature_;
            }

            //- Lists of connected featureEdges. Indices into featureEdges.
            const labelListList& featureSegments() const
            {
                return featureSegments_;
            }

            //- Indices into edges of additional edges.
            const labelList& extraEdges() const
            {
                return extraEdges_;
            }


        // Edit

            //- Read from boundaryMesh of polyMesh.
            void read(const polyMesh&);

            //- Read from triSurface
            void readTriSurface(const fileName&);

            //- Write to file.
            void writeTriSurface(const fileName&) const;

            //- Get bMesh index of nearest face for every boundary face in
            //  pMesh. Gets passed initial search box. If not found
            //  returns -1 for the face.
            labelList getNearest
            (
                const primitiveMesh& pMesh,
                const vector& searchSpan
            ) const;

            //- Take over patches onto polyMesh from nearest face in *this
            //  (from call to getNearest). Insert as
            //      -new set of patches (newMesh.addPatches)
            //      -topoChanges to change faces.
            // nearest is list of nearest face in *this for every boundary
            // face. oldPatches is list of existing patches in mesh.
            // newMesh is the mesh to which the new patches are added.
            // (so has to be constructed without patches).
            void patchify
            (
                const labelList& nearest,
                const polyBoundaryMesh& oldPatches,
                polyMesh& newMesh
            ) const;

        // Patches

            //- Get index of patch face is in
            label whichPatch(const label faceI) const;

            //- Get index of patch by name
            label findPatchID(const word& patchName) const;

            //- Get names of patches
            wordList patchNames() const;

            //- Add to back of patch list.
            void addPatch(const word& patchName);

            //- Delete from patch list.
            void deletePatch(const word& patchName);

            //- Change patch.
            void changePatchType(const word& patchName, const word& type);

            //- Recalculate face ordering and patches. Return old to new
            //  mapping.
            void changeFaces(const labelList& patchIDs, labelList& oldToNew);


        // Edges

            //- Set featureEdges, edgeToFeature, featureSegments according
            //  to angle of faces across edge
            void setFeatureEdges(const scalar minCos);

            //- Set extraEdges to edges 'near' to edgeI. Uses point-edge walk
            //  to determine 'near'.
            void setExtraEdges(const label edgeI);


        // Faces

            //- Simple triangulation of face subset. Returns number of triangles
            //  needed.
            label getNTris(const label faceI) const;

            //- Simple triangulation of face subset. TotalNTris is total number
            //  of triangles, nTris is per face number of triangles.
            label getNTris
            (
                const label startFaceI,
                const label nFaces,
                labelList& nTris
            ) const;

            //- Simple triangulation of face subset. TotalNTris is total number
            //  of triangles (from call to getNTris)
            //  triVerts is triangle vertices, three per triangle.
            void triangulate
            (
                const label startFaceI,
                const label nFaces,
                const label totalNTris,
                labelList& triVerts
            ) const;

            //- Number of points used in face subset.
            label getNPoints(const label startFaceI, const label nFaces) const;

            //- Same as triangulate but in local vertex numbering.
            //  (Map returned).
            void triangulateLocal
            (
                const label startFaceI,
                const label nFaces,
                const label totalNTris,
                labelList& triVerts,
                labelList& localToGlobal
            ) const;

        // Other

            // Flood filling without crossing protected edges.
            void markFaces
            (
                const labelList& protectedEdges,
                const label faceI,
                boolList& visited
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
