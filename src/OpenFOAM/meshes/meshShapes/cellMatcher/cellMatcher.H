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
    Foam::cellMatcher

Description
    Base class for cellshape matchers (hexMatch, prismMatch, etc.). These are
    classes which given a mesh and cell number find out the orientation of
    the cellShape and construct cell-vertex to mesh-vertex mapping and
    cell-face to mesh-face mapping.

    For example,
    \verbatim
        hexMatcher hex(mesh);
        cellShape shape;
        ..
        bool isHex = hex.match(celli, shape);
    \endverbatim
    Now shape is set to the correct Hex cellShape (if \a isHex is true)

    Alternatively there is direct access to the vertex and face mapping:
    \verbatim
        const labelList& hexVertLabels = hex.vertLabels();
        const labelList& hexFaceLabels = hex.faceLabels();
    \endverbatim
    Now
      - \c hexVertLabels[n] is vertex label of hex vertex n
      - \c hexFaceLabels[n] is face   label of hex vertex n

    Process of cellShape recognition consists of following steps:
    - renumber vertices of cell to local vertex numbers
    - construct (local to cell) addressing edge-to-faces
    - construct (local to cell) addressing vertex and face to index in face
    - find most unique face shape (e.g. triangle for prism)
    - walk (following either vertices in face or jumping from face to other
      face) to other faces and checking face sizes.
    - if necessary try other rotations of this face
      (only necessary for wedge, tet-wedge)
    - if necessary try other faces which most unique face shape
      (never necessary for hex degenerates)

    The whole calculation is done such that no lists are allocated during
    cell checking. E.g. localFaces_ are always sized to hold max. number
    of possible face vertices and a separate list is filled which holds
    the actusl face sizes.

    For now all hex-degenerates implemented. Numbering taken from picture in
    demoGuide.

SourceFiles
    cellMatcherI.H
    cellMatcher.C

\*---------------------------------------------------------------------------*/

#ifndef cellMatcher_H
#define cellMatcher_H

#include "labelList.H"
#include "faceList.H"
#include "boolList.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class primitiveMesh;
class cell;
class cellShape;
class cellModel;

/*---------------------------------------------------------------------------*\
                         Class cellMatcher Declaration
\*---------------------------------------------------------------------------*/

class cellMatcher
{
protected:

    // Static functions

        //- Given start and end of edge generate unique key
        inline static label edgeKey
        (
            const label numVert,
            const label v0,
            const label v1
        );

        //- Step along face either in righthand or lefthand direction
        inline static label nextVert(const label, const label, const bool);


    // Protected data

        // Map from mesh to local vertex numbering
        Map<label> localPoint_;

        //- Faces using local vertex numbering
        faceList localFaces_;

        //- Number of vertices per face in localFaces_
        labelList faceSize_;

        //- Map from local to mesh vertex numbering
        labelList pointMap_;

        //- Map from local to mesh face numbering
        labelList faceMap_;

        //- Map from 'edge' to neighbouring faces
        labelList edgeFaces_;

        //- pointFaceIndex[localVertI][localFacei] is index in localFace
        //  where localVertI is.
        labelListList pointFaceIndex_;

        //- After matching: holds mesh vertices in cellmodel order
        labelList vertLabels_;

        //- After matching: holds mesh faces in cellmodel order
        labelList faceLabels_;

        //- CellModel name
        const word cellModelName_;

        mutable const cellModel* cellModelPtr_;


    // Protected Member Functions

        //- Calculates localFaces. Returns number of local vertices (or -1
        //  if more than vertPerCell).
        label calcLocalFaces(const faceList& faces, const labelList& myFaces);

        //- Fill edge (start, end) to face number
        void calcEdgeAddressing(const label numVert);

        //- Fill vertex/face to index in face data structure
        void calcPointFaceIndex();

        //- Given start,end of edge lookup both faces sharing it and return
        //  face != localFacei
        label otherFace
        (
            const label numVert,
            const label v0,
            const label v1,
            const label localFacei
        ) const;


public:

    // Constructors

        //- Construct given mesh and shape factors
        cellMatcher
        (
            const label vertPerCell,
            const label facePerCell,
            const label maxVertPerFace,
            const word& cellModelName
        );

        //- Disallow default bitwise copy construction
        cellMatcher(const cellMatcher&) = delete;


    //- Destructor
    virtual ~cellMatcher()
    {}


    // Member Functions

        // Access

            inline const Map<label>& localPoint() const;
            inline const faceList& localFaces() const;
            inline const labelList& faceSize() const;
            inline const labelList& pointMap() const;
            inline const labelList& faceMap() const;
            inline const labelList& edgeFaces() const;
            inline const labelListList& pointFaceIndex() const;
            inline const labelList& vertLabels() const;
            inline const labelList& faceLabels() const;
            inline const cellModel& model() const;


        // Write

            void write(Ostream& os) const;

        // Cell shape dependent

            virtual label nVertPerCell() const = 0;

            virtual label nFacePerCell() const = 0;

            virtual label nMaxVertPerFace() const = 0;

            //- Hash value of all face sizes of this shape. Can be used for
            //  quick initial recognition.
            virtual label faceHashValue() const = 0;

            //- Check whether number of face sizes match the shape.
            virtual bool faceSizeMatch(const faceList&, const labelList&)
                const = 0;

            //- Low level shape recognition. Return true if matches.
            //  Works in detection mode only (checkOnly=true) or in exact
            //  matching. Returns true and sets vertLabels_.
            //  Needs faces, faceOwner of all faces in 'mesh' and cell number
            //  and labels of faces for this cell.
            //  celli only used in combination with faceOwner to detect owner
            //  status.
            virtual bool matchShape
            (
                const bool checkOnly,
                const faceList& faces,
                const labelList& faceOwner,
                const label celli,
                const labelList& myFaces
            ) = 0;

            //- Exact match. Uses faceSizeMatch.
            //  Returns true if cell matches shape exactly.
            virtual bool isA(const primitiveMesh& mesh, const label celli) = 0;

            //- Exact match given all the faces forming a cell. No checks
            //  on whether faces match up and form a closed shape.
            virtual bool isA(const faceList&) = 0;

            //- Like isA but also constructs a cellShape (if shape matches)
            virtual bool matches
            (
                const primitiveMesh& mesh,
                const label celli,
                cellShape& shape
            ) = 0;


    // Member Operators

        void operator=(const cellMatcher&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellMatcherI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
