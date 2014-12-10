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
    Foam::mapDistributePolyMesh

Description
    Class containing mesh-to-mesh mapping information after a mesh distribution
    where we send parts of meshes (using subsetting) to other processors
    and receive and reconstruct mesh.

    We store mapping from the bits-to-send to the complete starting mesh
    (subXXXMap) and from the received bits to their location in the new
    mesh (constructXXXMap).

SourceFiles
    mapDistributePolyMesh.C

\*---------------------------------------------------------------------------*/

#ifndef mapDistributePolyMesh_H
#define mapDistributePolyMesh_H

#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class mapDistributePolyMesh Declaration
\*---------------------------------------------------------------------------*/

class mapDistributePolyMesh
{
    // Private data

        const polyMesh& mesh_;

        //- Number of old live points
        const label nOldPoints_;

        //- Number of old live faces
        const label nOldFaces_;

        //- Number of old live cells
        const label nOldCells_;

        //- List of the old patch sizes
        labelList oldPatchSizes_;

        //- List of the old patch start labels
        const labelList oldPatchStarts_;

        //- List of numbers of mesh points per old patch
        const labelList oldPatchNMeshPoints_;


        //- Point distribute map
        const mapDistribute pointMap_;

        //- Face distribute map
        const mapDistribute faceMap_;

        //- Cell distribute map
        const mapDistribute cellMap_;

        //- Patch distribute map
        const mapDistribute patchMap_;



    // Private Member Functions

        void calcPatchSizes();

        //- Disallow default bitwise copy construct
        mapDistributePolyMesh(const mapDistributePolyMesh&);

        //- Disallow default bitwise assignment
        void operator=(const mapDistributePolyMesh&);


public:

    // Constructors

        //- Construct from components. Note that mesh has to be changed already
        //  since uses mesh.nPoints etc as the new size.
        mapDistributePolyMesh
        (
            const polyMesh& mesh,

            // mesh before changes
            const label nOldPoints,
            const label nOldFaces,
            const label nOldCells,
            const Xfer<labelList>& oldPatchStarts,
            const Xfer<labelList>& oldPatchNMeshPoints,

            // how to subset pieces of mesh to send across
            const Xfer<labelListList>& subPointMap,
            const Xfer<labelListList>& subFaceMap,
            const Xfer<labelListList>& subCellMap,
            const Xfer<labelListList>& subPatchMap,

            // how to reconstruct received mesh
            const Xfer<labelListList>& constructPointMap,
            const Xfer<labelListList>& constructFaceMap,
            const Xfer<labelListList>& constructCellMap,
            const Xfer<labelListList>& constructPatchMap
        );


    // Member Functions

        // Access

            const polyMesh& mesh() const
            {
                return mesh_;
            }

            //- Number of points in mesh before distribution
            label nOldPoints() const
            {
                return nOldPoints_;
            }

            //- Number of faces in mesh before distribution
            label nOldFaces() const
            {
                return nOldFaces_;
            }

            //- Number of cells in mesh before distribution
            label nOldCells() const
            {
                return nOldCells_;
            }

            //- List of the old patch sizes
            const labelList& oldPatchSizes() const
            {
                return oldPatchSizes_;
            }

            //- List of the old patch start labels
            const labelList& oldPatchStarts() const
            {
                return oldPatchStarts_;
            }

            //- List of numbers of mesh points per old patch
            const labelList& oldPatchNMeshPoints() const
            {
                return oldPatchNMeshPoints_;
            }

            //- Point distribute map
            const mapDistribute& pointMap() const
            {
                return pointMap_;
            }

            //- Face distribute map
            const mapDistribute& faceMap() const
            {
                return faceMap_;
            }

            //- Cell distribute map
            const mapDistribute& cellMap() const
            {
                return cellMap_;
            }

            //- Patch distribute map
            const mapDistribute& patchMap() const
            {
                return patchMap_;
            }


        // Edit

            //- distribute list of point data
            template<class T>
            void distributePointData(List<T>& lst) const
            {
                pointMap_.distribute(lst);
            }

            //- distribute list of face data
            template<class T>
            void distributeFaceData(List<T>& lst) const
            {
                faceMap_.distribute(lst);
            }

            //- distribute list of cell data
            template<class T>
            void distributeCellData(List<T>& lst) const
            {
                cellMap_.distribute(lst);
            }

            //- distribute list of patch data
            template<class T>
            void distributePatchData(List<T>& lst) const
            {
                patchMap_.distribute(lst);
            }


            //- distribute list of point/face/cell/patch indices.
            //  (Converts to boolList, distributes boolList and reconstructs)
            void distributePointIndices(labelList& pointIDs) const;

            void distributeFaceIndices(labelList& faceIDs) const;
            void distributeCellIndices(labelList& cellIDs) const;
            void distributePatchIndices(labelList& patchIDs) const;


            //- Correct for topo change.
            void updateMesh(const mapPolyMesh&)
            {
                notImplemented
                (
                    "mapDistributePolyMesh::updateMesh(const mapPolyMesh&)"
                );
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
