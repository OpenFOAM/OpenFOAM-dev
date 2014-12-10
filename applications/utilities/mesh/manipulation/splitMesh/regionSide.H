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
    Foam::regionSide

Description
    Determines the 'side' for every face and connected to a
    singly-connected (through edges) region of faces. Gets set of faces and
    a list of mesh edges ('fenceEdges') which should not be crossed.
    Used in splitting a mesh region.

    Determines:
    - For every face on the surface: whether the owner was visited
      from starting face.
    - List of faces using an internal point of the region visitable by
      edge-face-edge walking from the correct side of the region.

SourceFiles
    regionSide.C

\*---------------------------------------------------------------------------*/

#ifndef regionSide_H
#define regionSide_H

#include "HashSet.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                           Class regionSide Declaration
\*---------------------------------------------------------------------------*/

class regionSide
{
    // Private data

        //- For every face on region tells whether the owner is on the
        //  'regionside'.
        labelHashSet sideOwner_;

        //- Contains the faces using an internal point and visited face
        labelHashSet insidePointFaces_;


    // Private Member Functions

        //- Step across point to other edge on face
        static label otherEdge
        (
            const primitiveMesh& mesh,
            const label faceI,
            const label edgeI,
            const label pointI
        );

        //- From faceI, side cellI, cross to other faces/cells by
        //  face-cell walking and store visited faces and update sideOwner_.
        void visitConnectedFaces
        (
            const primitiveMesh& mesh,
            const labelHashSet& region,
            const labelHashSet& fenceEdges,
            const label cellI,
            const label faceI,
            labelHashSet& visitedFace
        );

        //- From edge on face connected to point on region (regionPointI) cross
        //  to all other edges using this point by walking across faces
        //  Does not cross regionEdges so stays on one side of region
        void walkPointConnectedFaces
        (
            const primitiveMesh& mesh,
            const labelHashSet& regionEdges,
            const label regionPointI,
            const label startFaceI,
            const label startEdgeI,
            labelHashSet& visitedEdges
        );

        //- Visits all internal points on region and marks edges reachable
        //  from sideOwner side (using walkPointConnectedFaces)
        void walkAllPointConnectedFaces
        (
            const primitiveMesh& mesh,
            const labelHashSet& regionFaces,
            const labelHashSet& fenceEdges
        );

public:

    //- Runtime type information
    ClassName("regionSide");

    // Static Functions

        //- Step across edge onto other face on cell
        static label otherFace
        (
            const primitiveMesh& mesh,
            const label cellI,
            const label excludeFaceI,
            const label edgeI
        );


    // Constructors

        //- Construct from components
        regionSide
        (
            const primitiveMesh& mesh,
            const labelHashSet& region,
            const labelHashSet& fenceEdges,   // labels of fence edges
            const label startCell,
            const label startFace
        );


    // Member Functions

        // Access

            const labelHashSet& sideOwner() const
            {
                return sideOwner_;
            }

            const labelHashSet& insidePointFaces() const
            {
                return insidePointFaces_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
