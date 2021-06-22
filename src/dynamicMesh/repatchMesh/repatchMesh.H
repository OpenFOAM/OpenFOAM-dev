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
    Foam::repatchMesh

Description
    Addressing for all faces on surface of mesh. Can either be read
    from polyMesh or from triSurface. Used for repatching existing meshes.

SourceFiles
    repatchMesh.C

\*---------------------------------------------------------------------------*/

#ifndef repatchMesh_H
#define repatchMesh_H

#include "repatchPatch.H"
#include "PrimitivePatch.H"
#include "PtrList.H"
#include "polyPatchList.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                        Class repatchMesh Declaration
\*---------------------------------------------------------------------------*/

class repatchMesh
{
public:

    // Public Typedefs

        //- Type of the mesh
        typedef PrimitivePatch<faceList, const pointField> rMesh;


private:

    // Static data

        //- Normal along which to divide faces into categories
        //  (used in getNearest)
        static const vector splitNormal_;

        //- Distance to face tolerance for getNearest. Triangles are considered
        //  near if they are nearer than distanceTol_*typDim where typDim is
        //  the largest distance from face centre to one of its vertices.
        static const scalar distanceTol_;


    // Private Data

        //- All boundary mesh data. Reconstructed every time faces are repatched
        autoPtr<rMesh> meshPtr_;

        //- Patches. Reconstructed every time faces are repatched.
        PtrList<repatchPatch> patches_;

        //- For every face in mesh() gives corresponding polyMesh face
        // (not sensible if mesh read from triSurface)
        labelList meshFace_;

        //- Points referenced by feature edges.
        pointField featurePoints_;

        //- Feature edges. Indices into featurePoints.
        edgeList featureEdges_;

        //- From feature edge to mesh edge.
        labelList featureToEdge_;

        //- From mesh edges to featureEdges_;
        labelList edgeToFeature_;

        //- Feature 'segments'. Collections of connected featureEdges.
        //  Indices into featureEdges_.
        labelListList featureSegments_;


    // Private Member Functions

        //- Number of connected feature edges.
        label nFeatureEdges(label pointi) const;

        //- Step to next feature edge
        label nextFeatureEdge(const label edgeI, const label vertI) const;

        //- Return connected list of feature edges.
        labelList collectSegment
        (
            const boolList& isFeaturePoint,
            const label startEdgeI,
            boolList& featVisited
        ) const;

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

        //- Finds area, starting at facei, delimited by borderEdge. Marks all
        //  faces thus visited with currentZone.
        void markZone
        (
            const boolList& borderEdge,
            label facei,
            label currentZone,
            labelList& faceZone
        ) const;


public:

    //- Runtime type information
    ClassName("repatchMesh");


    // Constructors

        //- Construct null
        repatchMesh();

        //- Disallow default bitwise copy construction
        repatchMesh(const repatchMesh&) = delete;


    //- Destructor
    ~repatchMesh();


    // Member Functions

        // Access

            //- Access the boundary mesh
            const rMesh& mesh() const
            {
                if (!meshPtr_.valid())
                {
                    FatalErrorInFunction
                        << "No mesh available. Probably mesh not yet"
                        << " read." << abort(FatalError);
                }
                return *meshPtr_;
            }

            //- Access the patches
            const PtrList<repatchPatch>& patches() const
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


        // Edit

            //- Read from repatchMesh of polyMesh.
            void read(const polyMesh&);

            //- Read from triSurface
            void readTriSurface(const fileName&);

            //- Get rMesh index of nearest face for every boundary face in
            //  pMesh. Gets passed initial search box. If not found
            //  returns -1 for the face.
            labelList getNearest
            (
                const primitiveMesh& pMesh,
                const vector& searchSpan
            ) const;


        // Patches

            //- Get index of patch face is in
            label whichPatch(const label facei) const;

            //- Get index of patch by name
            label findPatchID(const word& patchName) const;

            //- Add to back of patch list.
            void addPatch(const word& patchName);

            //- Delete from patch list.
            void deletePatch(const word& patchName);

            //- Change patch.
            void changePatchType(const word& patchName, const word& type);


        // Edges

            //- Set featureEdges, edgeToFeature, featureSegments according
            //  to angle of faces across edge
            void setFeatureEdges(const scalar minCos);


        // Other

            // Flood filling without crossing protected edges.
            void markFaces
            (
                const labelList& protectedEdges,
                const label facei,
                boolList& visited
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const repatchMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
