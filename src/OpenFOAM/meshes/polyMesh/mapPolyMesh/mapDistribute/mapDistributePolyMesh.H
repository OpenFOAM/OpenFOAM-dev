/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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


// Forward declaration of friend functions and operators

class mapDistributePolyMesh;

Istream& operator>>(Istream&, mapDistributePolyMesh&);
Ostream& operator<<(Ostream&, const mapDistributePolyMesh&);


/*---------------------------------------------------------------------------*\
                    Class mapDistributePolyMesh Declaration
\*---------------------------------------------------------------------------*/

class mapDistributePolyMesh
{
    // Private Data

        //- Number of old live points
        label nOldPoints_;

        //- Number of old live faces
        label nOldFaces_;

        //- Number of old live cells
        label nOldCells_;

        //- List of the old patch sizes
        labelList oldPatchSizes_;

        //- List of the old patch start labels
        labelList oldPatchStarts_;

        //- List of numbers of mesh points per old patch
        labelList oldPatchNMeshPoints_;


        //- Point distribute map
        mapDistribute pointMap_;

        //- Face distribute map
        mapDistribute faceMap_;

        //- Cell distribute map
        mapDistribute cellMap_;

        //- Patch distribute map
        mapDistribute patchMap_;


    // Private Member Functions

        void calcPatchSizes();


public:

    // Constructors

        //- Construct null
        mapDistributePolyMesh();

        //- Move constructor from components.
        //  Note that mesh has to be changed already
        //  since uses mesh.nPoints etc as the new size.
        mapDistributePolyMesh
        (
            const polyMesh& mesh,

            // mesh before changes
            const label nOldPoints,
            const label nOldFaces,
            const label nOldCells,
            labelList&& oldPatchStarts,
            labelList&& oldPatchNMeshPoints,

            // how to subset pieces of mesh to send across
            labelListList&& subPointMap,
            labelListList&& subFaceMap,
            labelListList&& subCellMap,
            labelListList&& subPatchMap,

            // how to reconstruct received mesh
            labelListList&& constructPointMap,
            labelListList&& constructFaceMap,
            labelListList&& constructCellMap,
            labelListList&& constructPatchMap,

            const bool subFaceHasFlip = false,
            const bool constructFaceHasFlip = false
        );

        //- Move constructor from components
        mapDistributePolyMesh
        (
            // mesh before changes
            const label nOldPoints,
            const label nOldFaces,
            const label nOldCells,
            labelList&& oldPatchStarts,
            labelList&& oldPatchNMeshPoints,

            // how to subset pieces of mesh to send across
            mapDistribute&& pointMap,
            mapDistribute&& faceMap,
            mapDistribute&& cellMap,
            mapDistribute&& patchMap
        );

        //- Move constructor
        mapDistributePolyMesh(mapDistributePolyMesh&&);

        //- Construct from Istream
        mapDistributePolyMesh(Istream&);

        //- Disallow default bitwise copy construction
        mapDistributePolyMesh(const mapDistributePolyMesh&) = delete;


    // Member Functions

        // Access

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


        // Other

            //- Transfer the contents of the argument and annul the argument.
            void transfer(mapDistributePolyMesh&);

            //- Distribute list of point data
            template<class T>
            void distributePointData(List<T>& lst) const
            {
                pointMap_.distribute(lst);
            }

            //- Distribute list of face data
            template<class T>
            void distributeFaceData(List<T>& lst) const
            {
                faceMap_.distribute(lst);
            }

            //- Distribute list of cell data
            template<class T>
            void distributeCellData(List<T>& lst) const
            {
                cellMap_.distribute(lst);
            }

            //- Distribute list of patch data
            template<class T>
            void distributePatchData(List<T>& lst) const
            {
                patchMap_.distribute(lst);
            }


            //- Distribute list of point/face/cell/patch indices.
            //  (Converts to boolList, distributes boolList and reconstructs)
            void distributePointIndices(labelList& pointIDs) const;

            void distributeFaceIndices(labelList& faceIDs) const;
            void distributeCellIndices(labelList& cellIDs) const;
            void distributePatchIndices(labelList& patchIDs) const;


            //- Correct for topo change.
            void updateMesh(const mapPolyMesh&)
            {
                NotImplemented;
            }


    // Member Operators

        void operator=(const mapDistributePolyMesh&);

        void operator=(mapDistributePolyMesh&&);


    // IOstream Operators

        //- Read dictionary from Istream
        friend Istream& operator>>(Istream&, mapDistributePolyMesh&);

        //- Write dictionary to Ostream
        friend Ostream& operator<<(Ostream&, const mapDistributePolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
