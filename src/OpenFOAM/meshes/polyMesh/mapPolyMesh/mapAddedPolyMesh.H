/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::mapAddedPolyMesh

Description
    Class containing mesh-to-mesh mapping information after a mesh addition
    where we add a mesh ('added mesh') to an old mesh, creating a new mesh.

    We store mapping from the old to the new mesh and from the added mesh
    to the new mesh.

    Note: Might need some more access functions or maybe some zone maps?

SourceFiles
    mapAddedPolyMesh.C

\*---------------------------------------------------------------------------*/

#ifndef mapAddedPolyMesh_H
#define mapAddedPolyMesh_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class mapAddedPolyMesh Declaration
\*---------------------------------------------------------------------------*/

class mapAddedPolyMesh
{
    // Private data

        //- Old mesh points/face/cells
        label nOldPoints_;
        label nOldFaces_;
        label nOldCells_;

        //- Added mesh points/faces/cells
        label nAddedPoints_;
        label nAddedFaces_;
        label nAddedCells_;


        //- From old mesh points to new points
        labelList oldPointMap_;
        //- From old mesh faces to new faces
        labelList oldFaceMap_;
        //- From old mesh cells to new cells
        labelList oldCellMap_;

        //- From added mesh points to new points
        labelList addedPointMap_;
        //- From added mesh faces to new faces
        labelList addedFaceMap_;
        //- From added mesh cells to new cells
        labelList addedCellMap_;

        //- Original mesh to new mesh patch map. -1 for deleted patches.
        labelList oldPatchMap_;

        //- Added mesh to new mesh patch map. -1 for deleted patches.
        labelList addedPatchMap_;

        //- Original patch sizes on old mesh
        labelList oldPatchSizes_;

        //- Original patch starts
        labelList oldPatchStarts_;


public:

    // Constructors

        //- Construct from components
        mapAddedPolyMesh
        (
            const label nOldPoints,
            const label nOldFaces,
            const label nOldCells,
            const label nAddedPoints,
            const label nAddedFaces,
            const label nAddedCells,
            const labelList& oldPointMap,
            const labelList& oldFaceMap,
            const labelList& oldCellMap,

            const labelList& addedPointMap,
            const labelList& addedFaceMap,
            const labelList& addedCellMap,

            const labelList& oldPatchMap,
            const labelList& addedPatchMap,
            const labelList& oldPatchSizes,
            const labelList& oldPatchStarts
        );


    // Member Functions

        // Access

            // Old mesh data

                label nOldPoints() const
                {
                    return nOldPoints_;
                }

                label nOldFaces() const
                {
                    return nOldFaces_;
                }

                label nOldCells() const
                {
                    return nOldCells_;
                }


                //- From old mesh point/face/cell to new mesh point/face/cell.
                const labelList& oldPointMap() const
                {
                    return oldPointMap_;
                }
                const labelList& oldFaceMap() const
                {
                    return oldFaceMap_;
                }
                const labelList& oldCellMap() const
                {
                    return oldCellMap_;
                }

                //- From old patch index to new patch index or -1 if patch
                //  not present (since 0 size)
                const labelList& oldPatchMap() const
                {
                    return oldPatchMap_;
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

                //- Number of old internal faces
                label nOldInternalFaces() const
                {
                    return oldPatchStarts_[0];
                }


            // Added mesh data

                label nAddedPoints() const
                {
                    return nAddedPoints_;
                }

                label nAddedFaces() const
                {
                    return nAddedFaces_;
                }

                label nAddedCells() const
                {
                    return nAddedCells_;
                }

                //- From added mesh point/face/cell to new mesh point/face/cell.
                const labelList& addedPointMap() const
                {
                    return addedPointMap_;
                }
                const labelList& addedFaceMap() const
                {
                    return addedFaceMap_;
                }
                const labelList& addedCellMap() const
                {
                    return addedCellMap_;
                }

                //- From added mesh patch index to new patch index or -1 if
                //  patch not present (since 0 size)
                const labelList& addedPatchMap() const
                {
                    return addedPatchMap_;
                }


        // Edit

            void updateMesh(const mapPolyMesh&)
            {
                NotImplemented;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
