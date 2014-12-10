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
    Foam::meshCutAndRemove

Description
    like meshCutter but also removes non-anchor side of cell.

SourceFiles
    meshCutAndRemove.C

\*---------------------------------------------------------------------------*/

#ifndef meshCutAndRemove_H
#define meshCutAndRemove_H

#include "edgeVertex.H"
#include "boolList.H"
#include "labelList.H"
#include "typeInfo.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
class polyTopoChange;
class cellCuts;
class polyMesh;
class face;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class meshCutAndRemove Declaration
\*---------------------------------------------------------------------------*/

class meshCutAndRemove
:
    public edgeVertex
{
    // Private data

        //- Faces added in last setRefinement. Per split cell label of added
        //  face
        Map<label> addedFaces_;

        //- Points added in last setRefinement. Per split edge label of added
        //  point
        HashTable<label, edge, Hash<edge> > addedPoints_;


    // Private Static Functions

        // Returns -1 or index in elems1 of first shared element.
        static label firstCommon(const labelList& lst1, const labelList& lst2);

        //- Do the elements of edge appear in consecutive order in the list
        static bool isIn(const edge&, const labelList&);


    // Private Member Functions

        //- Returns -1 or the cell in cellLabels that is cut.
        label findCutCell(const cellCuts&, const labelList&) const;

        //- Returns first pointI in pointLabels that uses an internal
        //  face. Used to find point to inflate cell/face from (has to be
        //  connected to internal face)
        label findInternalFacePoint(const labelList& pointLabels) const;

        //- Find point on face that is part of original mesh and that is
        //  point connected to the patch
        label findPatchFacePoint(const face& f, const label patchI) const;

        //- Get new owner and neighbour of face. Checks anchor points to see if
        // need to get original or added cell.
        void faceCells
        (
            const cellCuts& cuts,
            const label exposedPatchI,
            const label faceI,
            label& own,
            label& nei,
            label& patchID
        ) const;

        //- Get zone information for face.
        void getZoneInfo
        (
            const label faceI,
            label& zoneID,
            bool& zoneFlip
        ) const;

        //- Adds a face from point. Flips face if owner>neighbour
        void addFace
        (
            polyTopoChange& meshMod,
            const label faceI,
            const label masterPointI,
            const face& newFace,
            const label owner,
            const label neighbour,
            const label patchID
        );

        //- Modifies existing faceI for either new owner/neighbour or
        //  new face points. Checks if anything changed and flips face
        //  if owner>neighbour
        void modFace
        (
            polyTopoChange& meshMod,
            const label faceI,
            const face& newFace,
            const label owner,
            const label neighbour,
            const label patchID
        );

        // Copies face starting from startFp. Jumps cuts. Marks visited
        // vertices in visited.
        void copyFace
        (
            const face& f,
            const label startFp,
            const label endFp,
            face& newFace
        ) const;

        //- Split face along cut into two faces. Faces are in same point
        //  order as original face (i.e. maintain normal direction)
        void splitFace
        (
            const face& f,
            const label v0,
            const label v1,

            face& f0,
            face& f1
        ) const;

        //- Add cuts of edges to face
        face addEdgeCutsToFace(const label faceI) const;

        //- Convert loop of cuts into face.
        face loopToFace
        (
            const label cellI,
            const labelList& loop
        ) const;



        //- Disallow default bitwise copy construct
        meshCutAndRemove(const meshCutAndRemove&);

        //- Disallow default bitwise assignment
        void operator=(const meshCutAndRemove&);

public:

    //- Runtime type information
    ClassName("meshCutAndRemove");


    // Constructors

        //- Construct from mesh
        meshCutAndRemove(const polyMesh& mesh);


    // Member Functions

        // Edit

            //- Do actual cutting with cut description. Inserts mesh changes
            //  into meshMod.
            //  cuts: all loops and topological information
            //  cutPatch: for every cell that has loop the patch number
            //  exposedPatch: patch for other exposed faces
            void setRefinement
            (
                const label exposedPatchI,
                const cellCuts& cuts,
                const labelList& cutPatch,
                polyTopoChange& meshMod
            );

            //- Force recalculation of locally stored data on topological change
            void updateMesh(const mapPolyMesh&);


        // Access

            //- Faces added. Per split cell label of added face
            const Map<label>& addedFaces() const
            {
                return addedFaces_;
            }

            //- Points added. Per split edge label of added point.
            //  (note: fairly useless across topology changes since one of the
            //  points of the edge will probably disappear)
            const HashTable<label, edge, Hash<edge> >& addedPoints() const
            {
                return addedPoints_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
