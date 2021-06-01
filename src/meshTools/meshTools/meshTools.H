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

Namespace
    Foam::meshTools

Description
    Collection of static functions to do various simple mesh related things.

SourceFiles
    meshTools.C

\*---------------------------------------------------------------------------*/

#ifndef meshTools_H
#define meshTools_H

#include "label.H"
#include "vector.H"
#include "triad.H"
#include "labelList.H"
#include "pointField.H"
#include "faceList.H"
#include "cellList.H"
#include "primitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class primitiveMesh;
class polyMesh;

/*---------------------------------------------------------------------------*\
                        Namespace meshTools Declaration
\*---------------------------------------------------------------------------*/

namespace meshTools
{
    // Bit identifiers for octants (p=plus, m=min e.g. plusXminYminZ)

    static const label mXmYmZ = 0;
    static const label pXmYmZ = 1;
    static const label mXpYmZ = 2;
    static const label pXpYmZ = 3;

    static const label mXmYpZ = 4;
    static const label pXmYpZ = 5;
    static const label mXpYpZ = 6;
    static const label pXpYpZ = 7;

    static const label mXmYmZMask = 1 << mXmYmZ;
    static const label pXmYmZMask = 1 << pXmYmZ;
    static const label mXpYmZMask = 1 << mXpYmZ;
    static const label pXpYmZMask = 1 << pXpYmZ;

    static const label mXmYpZMask = 1 << mXmYpZ;
    static const label pXmYpZMask = 1 << pXmYpZ;
    static const label mXpYpZMask = 1 << mXpYpZ;
    static const label pXpYpZMask = 1 << pXpYpZ;


    // Normal handling

        //- Check if n is in same direction as normals of all faceLabels
        bool visNormal
        (
            const vector& n,
            const vectorField& faceNormals,
            const labelList& faceLabels
        );

        //- Calculate point normals on a 'box' mesh (all edges aligned with
        //  coordinate axes)
        vectorField calcBoxPointNormals(const primitivePatch& pp);

        //- Normalised edge vector
        vector normEdgeVec(const primitiveMesh&, const label edgeI);


    // OBJ writing

        //- Write obj representation of point
        void writeOBJ
        (
            Ostream& os,
            const point& pt
        );

        //- Write obj representation of a triad. Requires the location of the
        //  triad to be supplied
        void writeOBJ
        (
            Ostream& os,
            const triad& t,
            const point& pt
        );

        //- Write obj representation of a line connecting two points
        //  Need to keep track of points that have been added. count starts at 0
        void writeOBJ
        (
            Ostream& os,
            const point& p1,
            const point& p2,
            label& count
        );

        //- Write obj representation of a point p1 with a vector from p1 to p2
        void writeOBJ
        (
            Ostream& os,
            const point& p1,
            const point& p2
        );

        //- Write obj representation of faces subset
        template<class FaceType>
        void writeOBJ
        (
            Ostream& os,
            const UList<FaceType>&,
            const pointField&,
            const labelList& faceLabels
        );

        //- Write obj representation of faces
        template<class FaceType>
        void writeOBJ
        (
            Ostream& os,
            const UList<FaceType>&,
            const pointField&
        );

        //- Write obj representation of cell subset
        void writeOBJ
        (
            Ostream& os,
            const cellList&,
            const faceList&,
            const pointField&,
            const labelList& cellLabels
        );


    // Cell/face/edge walking

        //- Is edge used by cell
        bool edgeOnCell
        (
            const primitiveMesh&,
            const label celli,
            const label edgeI
        );

        //- Is edge used by face
        bool edgeOnFace
        (
            const primitiveMesh&,
            const label facei,
            const label edgeI
        );

        //- Is face used by cell
        bool faceOnCell
        (
            const primitiveMesh&,
            const label celli,
            const label facei
        );

        //- Return edge among candidates that uses the two vertices.
        label findEdge
        (
            const edgeList& edges,
            const labelList& candidates,
            const label v0,
            const label v1
        );

        //- Return edge between two vertices. Returns -1 if no edge.
        label findEdge
        (
            const primitiveMesh&,
            const label v0,
            const label v1
        );

        //- Return edge shared by two faces. Throws error if no edge found.
        label getSharedEdge
        (
            const primitiveMesh&,
            const label f0,
            const label f1
        );

        //- Return face shared by two cells. Throws error if none found.
        label getSharedFace
        (
            const primitiveMesh&,
            const label cell0,
            const label cell1
        );

        //- Get faces on cell using edgeI. Throws error if no two found.
        void getEdgeFaces
        (
            const primitiveMesh&,
            const label celli,
            const label edgeI,
            label& face0,
            label& face1
        );

        //- Return label of other edge (out of candidates edgeLabels)
        //  connected to vertex but not edgeI. Throws error if none found.
        label otherEdge
        (
            const primitiveMesh&,
            const labelList& edgeLabels,
            const label edgeI,
            const label vertI
        );

        //- Return face on cell using edgeI but not facei. Throws error
        //  if none found.
        label otherFace
        (
            const primitiveMesh&,
            const label celli,
            const label facei,
            const label edgeI
        );

        //- Return cell on other side of face. Throws error
        //  if face not internal.
        label otherCell
        (
            const primitiveMesh&,
            const label celli,
            const label facei
        );

        //- Returns label of edge nEdges away from startEdge (in the direction
        // of startVertI)
        label walkFace
        (
            const primitiveMesh&,
            const label facei,
            const label startEdgeI,
            const label startVertI,
            const label nEdges
        );


    // Constraints on position

        //- Set the constrained components of position to mesh centre
        void constrainToMeshCentre
        (
            const polyMesh& mesh,
            point& pt
        );
        void constrainToMeshCentre
        (
            const polyMesh& mesh,
            pointField& pt
        );

        //- Set the constrained components of directions/velocity to zero
        void constrainDirection
        (
            const polyMesh& mesh,
            const Vector<label>& dirs,
            vector& d
        );
        void constrainDirection
        (
            const polyMesh& mesh,
            const Vector<label>& dirs,
            vectorField& d
        );


    // Hex only functionality.

        //- Given edge on hex find other 'parallel', non-connected edges.
        void getParallelEdges
        (
            const primitiveMesh&,
            const label celli,
            const label e0,
            label&,
            label&,
            label&
        );

        //- Given edge on hex find all 'parallel' (i.e. non-connected)
        //  edges and average direction of them
        vector edgeToCutDir
        (
            const primitiveMesh&,
            const label celli,
            const label edgeI
        );

        //- Reverse of edgeToCutDir: given direction find edge bundle and
        //  return one of them.
        label cutDirToEdge
        (
            const primitiveMesh&,
            const label celli,
            const vector& cutDir
        );

} // End namespace meshTools


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
