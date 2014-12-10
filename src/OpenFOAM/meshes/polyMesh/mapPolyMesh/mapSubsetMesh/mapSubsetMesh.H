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
    Foam::mapSubsetMesh

Description
    Class containing mesh-to-mesh mapping information after a subset operation

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef mapSubsetMesh_H
#define mapSubsetMesh_H

#include "refCount.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mapSubsetMesh Declaration
\*---------------------------------------------------------------------------*/

class mapSubsetMesh
:
    public refCount
{
    // Private data

        //- Number of old live points
        label nOldPoints_;

        //- Number of old live faces
        label nOldFaces_;

        //- Number of old live cells
        label nOldCells_;

        //- Point mapping array
        const labelList pointMap_;

        //- Face mapping array
        const labelList faceMap_;

        //- Cell mapping array
        const labelList cellMap_;

        //- Patch for exposed faces
        const label exposedPatchID_;

        //- List of the old patch sizes
        labelList oldPatchSizes_;

        //- List of the old patch start labels
        const labelList oldPatchStarts_;

        //- List of numbers of mesh points per old patch
        const labelList oldPatchNMeshPoints_;

public:

    // Constructors

        //- Construct from components
        mapSubsetMesh
        (
            const label nOldPoints,
            const label nOldFaces,
            const label nOldCells,
            const labelList& pointMap,
            const labelList& faceMap,
            const labelList& cellMap,
            const label exposedPatchID,
            const labelList& oldPatchStarts,
            const labelList& oldPatchNMeshPoints
        );

        //- Construct from components and optionally reuse storage
        mapSubsetMesh
        (
            const label nOldPoints,
            const label nOldFaces,
            const label nOldCells,
            labelList& pointMap,
            labelList& faceMap,
            labelList& cellMap,
            const label exposedPatchID,
            labelList& oldPatchStarts,
            labelList& oldPatchNMeshPoints,
            const bool reUse
        );


    // Member Functions

        // Access

            //- Number of old points
            label nOldPoints() const
            {
                return nOldPoints_;
            }

            //- Number of old internal faces
            label nOldInternalFaces() const
            {
                return oldPatchStarts_[0];
            }

            //- Number of old faces
            label nOldFaces() const
            {
                return nOldFaces_;
            }

            //- Number of old cells
            label nOldCells() const
            {
                return nOldCells_;
            }

            //- Patch that exposed faces were put into
            label exposedPatchID() const
            {
                return exposedPatchID_;
            }

            //- Old point map.
            //  Contains the old point label for all points of the subsetted
            //  mesh
            const labelList& pointMap() const
            {
                return pointMap_;
            }

            //- Old face map.
            //  Contains the old point label for all faces of the subsetted
            //  mesh
            const labelList& faceMap() const
            {
                return faceMap_;
            }

            //- Old cell map.
            //  Contains the old point label for all cells of the subsetted
            //  mesh
            const labelList& cellMap() const
            {
                return cellMap_;
            }

            //- Return list of the old patch sizes
            const labelList& oldPatchSizes() const
            {
                return oldPatchSizes_;
            }

            //- Return list of the old patch start labels
            const labelList& oldPatchStarts() const
            {
                return oldPatchStarts_;
            }

            //- Return numbers of mesh points per old patch
            const labelList& oldPatchNMeshPoints() const
            {
                return oldPatchNMeshPoints_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
