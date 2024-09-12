/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::domainDecomposition

Description
    Automatic domain decomposition class for finite-volume meshes

SourceFiles
    domainDecomposition.C
    domainDecompositionDecompose.C
    domainDecompositionReconstruct.C
    domainDecompositionReconstruct.C

\*---------------------------------------------------------------------------*/

#ifndef domainDecomposition_H
#define domainDecomposition_H

#include "fvMesh.H"
#include "processorRunTimes.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "regionName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class faceCoupleInfo;
class multiDomainDecomposition;

/*---------------------------------------------------------------------------*\
                     Class domainDecomposition Declaration
\*---------------------------------------------------------------------------*/

class domainDecomposition
{
    // Private Typedefs

        //- Table to map from a labelPair to a label
        typedef
            HashTable<label, labelPair, Hash<labelPair>>
            labelPairLabelTable;


    // Private Data

        //- Optimisation switch for reconstruction algorithm. See corresponding
        //  method below.
        static bool sortReconstructNonConformalCyclicAddressing_;

        //- Run times
        const processorRunTimes& runTimes_;

        //- Mesh path
        const fileName meshPath_;

        //- Region name
        const word regionName_;

        //- The complete mesh
        autoPtr<fvMesh> completeMesh_;

        //- The processor meshes
        PtrList<fvMesh> procMeshes_;

        //- The region meshes
        const multiDomainDecomposition& regionMeshes_;


        // Complete mesh to processor mesh addressing

            //- For each complete cell, the processor index
            labelList cellProc_;


        // Processor mesh to complete mesh addressing

            //- For each processor point, the complete point index
            labelListList procPointAddressing_;

            //- For each processor face, the complete face index
            // Note: Face turning index is stored as the sign on addressing
            // Only the processor boundary faces are affected: if the sign of
            // the index is negative, the processor face is the reverse of the
            // original face. In order to do this properly, all face
            // indices will be incremented by 1 and the decremented as
            // necessary to avoid the problem of face number zero having no
            // sign.
            labelListList procFaceAddressing_;

            //- For each processor cell, the complete point index
            labelListList procCellAddressing_;


        // Finite volume specific processor to complete mesh addressing

            //- Labels of finite volume faces for each processor boundary
            //  (constructed on demand)
            mutable PtrList<surfaceLabelField::Boundary> procFaceAddressingBf_;


    // Private Member Functions

        // Decomposition

            //- Mark all elements with value or -2 if occur twice
            static void mark
            (
                const labelList& zoneElems,
                const label zoneI,
                labelList& elementToZone
            );

            //- Add face to inter-processor patch
            void addInterProcFace
            (
                const label facei,
                const label ownerProc,
                const label nbrProc,
                const label subPatchID,
                List<Map<label>>&,
                List<DynamicList<DynamicList<label>>>&,
                List<labelListList>&,
                List<labelListList>&
            ) const;

            //- Call the decomposition method and return the processor index
            //  that each cell is being distributed to
            labelList distributeCells();

            //- Generate sub patch info for processor cyclics
            void processInterCyclics
            (
                const labelList& cellProc,
                const polyBoundaryMesh& patches,
                List<DynamicList<DynamicList<label>>>& interPatchFaces,
                List<Map<label>>& procNbrToInterPatch,
                List<labelListList>& subPatchIDs,
                List<labelListList>& subPatchStarts,
                bool owner,
                bool first
            ) const;

            //- Decompose the complete mesh to create the processor meshes and
            //  populate the addressing
            void decompose();

            //- Decompose the complete mesh's points and apply the result to
            //  the processor meshes
            void decomposePoints();


        // Reconstruction

            //- Find shared points and faces between two meshes that are to be
            //  added together
            static autoPtr<faceCoupleInfo> determineCoupledFaces
            (
                const label masterMeshProcStart,
                const label masterMeshProcEnd,
                const polyMesh& masterMesh,
                const label meshToAddProcStart,
                const label meshToAddProcEnd,
                const polyMesh& meshToAdd
            );

            //- Reconstruct the processor meshes to create the complete mesh
            //  and populate the addressing
            void reconstruct();

            //- Reconstruct the processor meshes' points and apply the result
            //  to the complete mesh
            void reconstructPoints();


        // Non conformal

            //- Is the complete mesh conformal?
            bool completeConformal() const;

            //- Are the processor meshes conformal?
            bool procsConformal() const;

            //- Generate the reverse of proc-face-face addressing. -1 indicates
            //  a face that is on a (conformal) processor boundary and hence
            //  has multiple associated proc-face indices.
            labelList completeFaceAddressing() const;

            //- Generate a map from non-conformal-cyclic patch index and an
            //  owner-neighbour pair of processor indices to the corresponding
            //  non-conformal processor-cyclic patch index
            List<labelPairLabelTable> nonConformalCyclicProcCyclics() const;

            //- Generate the processor offsets for the non-conformal mapped
            //  wall patches. Optionally extend the lists to nProcs+1 and add
            //  the total size at the end.
            PtrList<labelListList>
                nonConformalMappedWallProcOffsets(const bool appendSize) const;

            //- Decompose the addressing of a non-conformal-cyclic patch.
            //  Helper for nonConformalProcFaceAddressingBf below.
            void decomposeNonConformalCyclicAddressing
            (
                const label nccPatchi,
                const List<labelPairLabelTable>& nonConformalCyclicProcCyclics,
                List<List<DynamicList<label>>>& nonConformalProcFaceAddressingBf
            ) const;

            //- Decompose the addressing of a non-conformal-mapped-wall patch.
            //  Helper for nonConformalProcFaceAddressingBf below.
            void decomposeNonConformalMappedWallAddressing
            (
                const label ncmwPatchi,
                PtrList<labelListList>& nonConformalMappedWallProcOffsets,
                List<List<DynamicList<label>>>& nonConformalProcFaceAddressingBf
            ) const;

            //- Decompose the addressing of a non-conformal-error patch.
            //  Helper for nonConformalProcFaceAddressingBf below.
            void decomposeNonConformalErrorAddressing
            (
                const label ncePatchi,
                List<List<DynamicList<label>>>& nonConformalProcFaceAddressingBf
            ) const;

            //- Reconstruct the addressing of a non-conformal-cyclic patch.
            //  Helper for nonConformalProcFaceAddressingBf below.
            void reconstructNonConformalCyclicAddressing
            (
                const label nccPatchi,
                const List<labelPairLabelTable>& nonConformalCyclicProcCyclics,
                List<List<DynamicList<label>>>& nonConformalProcFaceAddressingBf
            ) const;

            //- Reconstruct the addressing of a non-conformal-cyclic patch.
            //  Helper for nonConformalProcFaceAddressingBf below.
            //
            //  !!! Alternative implementation based on a sort. Has complexity
            //  nFaces*log2(nFaces), whether as the "standard" function has
            //  complexity nProcs*nProcs*nFaces. It was thought that this
            //  version would be better, but testing hasn't found any cases for
            //  which this is the case. Enabled with an optimisation switch
            //  with the same name as the method.
            //
            void sortReconstructNonConformalCyclicAddressing
            (
                const label nccPatchi,
                const List<labelPairLabelTable>& nonConformalCyclicProcCyclics,
                List<List<DynamicList<label>>>& nonConformalProcFaceAddressingBf
            ) const;

            //- Reconstruct the addressing of a non-conformal-mapped-wall patch.
            //  Helper for nonConformalProcFaceAddressingBf below.
            void reconstructNonConformalMappedWallAddressing
            (
                const label ncmwPatchi,
                List<List<DynamicList<label>>>& nonConformalProcFaceAddressingBf
            ) const;

            //- Reconstruct the addressing of a non-conformal-error patch.
            //  Helper for nonConformalProcFaceAddressingBf below.
            void reconstructNonConformalErrorAddressing
            (
                const label ncePatchi,
                List<List<DynamicList<label>>>& nonConformalProcFaceAddressingBf
            ) const;

            //- Generate the procFaceAddressingBf for the non-conformal patches
            //  only. Helper for procFaceAddressingBf.
            List<List<DynamicList<label>>>
                nonConformalProcFaceAddressingBf() const;

            //- Map the polyFacesBf from the processor to the complete meshes
            //  and unconform the complete mesh
            void unconformComplete();

            //- Map the polyFacesBf from the complete to the processor meshes
            //  and unconform the processor meshes
            void unconformProcs();

            //- If the complete mesh is non-conformal and the processor meshes
            //  are not, then unconform the processor meshes. And vice versa.
            void unconform();


        //- Compare two instances. A return value of -1 means a is newer than b
        //  (i.e., the arguments are in reverse order), 0 means a is equal to
        //  b, and +1 means a is older than b (in order).
        label compareInstances(const fileName& a, const fileName& b) const;

        //- Validate that the complete mesh has been generated or read
        void validateComplete() const;

        //- Validate that the processor meshes have been generated or read
        void validateProcs() const;

        //- Read the complete mesh
        void readComplete();

        //- Read the processor meshes
        void readProcs();

        //- Read the addressing
        void readCompleteAddressing();

        //- Read the addressing
        void readProcsAddressing();

        //- Read the addressing
        void readAddressing();

        //- Read-update the complete and processor meshes for a change in time
        fvMesh::readUpdateState readUpdate();

        //- Write the decomposition addressing
        void writeCompleteAddressing() const;

        //- Write the decomposition addressing
        void writeProcsAddressing() const;

        //- Write the decomposition addressing
        void writeAddressing() const;

        //- Write the processor meshes' points for an instance different from
        //  the current. Used to additionally write out points associated with
        //  the face instance.
        void writeProcPoints(const fileName& inst);

        //- Write the complete mesh's points for an instance different from
        //  the current. Used to additionally write out points associated with
        //  the face instance.
        void writeCompletePoints(const fileName& inst);


public:

    //- Runtime type information
    TypeName("domainDecomposition");


    // Constructors

        //- Construct from processor run times and region name
        domainDecomposition
        (
            const processorRunTimes& procRunTimes,
            const fileName& meshPath,
            const word& regionName,
            const multiDomainDecomposition& regionMeshes =
                NullObjectRef<multiDomainDecomposition>()
        );


    //- Destructor
    virtual ~domainDecomposition();


    // Member Functions

        //- Access the run times
        inline const processorRunTimes& procRunTimes() const
        {
            return runTimes_;
        }

        //- Access the mesh path
        inline const fileName& meshPath() const
        {
            return meshPath_;
        }

        //- Access the region name
        inline const word& regionName() const
        {
            return regionName_;
        }

        //- Access the global mesh
        inline const fvMesh& completeMesh() const
        {
            validateComplete();
            return completeMesh_();
        }

        //- Access the processor meshes
        inline const PtrList<fvMesh>& procMeshes() const
        {
            validateProcs();
            return procMeshes_;
        }

        //- Return the number of processors in the decomposition
        inline label nProcs() const
        {
            return runTimes_.nProcs();
        }

        //- Read in the complete mesh. Read the processor meshes if they are
        //  available and up to date relative to the complete mesh, or
        //  decompose if not. Return whether or not decomposition happened.
        bool readDecompose();

        //- Post-read-construction steps for the meshes after read-decompose
        void postReadDecompose();

        //- Synchronise non-conformal structure after read-decompose
        void unconformReadDecompose();

        //- Write following read-decompose
        void writeReadDecompose(const bool decomposed, const bool doSets);

        //- Read in the processor meshes. Read the complete mesh if it is
        //  available and up to date relative to the processor meshes, or
        //  reconstruct if not. Return whether or not reconstruction happened.
        bool readReconstruct();

        //- Post-read-construction steps for the meshes after read-reconstruct
        void postReadReconstruct();

        //- Synchronise non-conformal structure after read-reconstruct
        void unconformReadReconstruct();

        //- Write following read-reconstruct
        void writeReadReconstruct(const bool reconstructed, const bool doSets);

        //- Read-update for decomposition
        fvMesh::readUpdateState readUpdateDecompose();

        //- Complete read-update-decompose
        void postReadUpdateDecompose(const fvMesh::readUpdateState stat);

        //- Synchronise non-conformal structure after read-update-decompose
        void unconformReadUpdateDecompose(const fvMesh::readUpdateState stat);

        //- Read-update for reconstruction
        fvMesh::readUpdateState readUpdateReconstruct();

        //- Complete read-update-reconstruct
        void postReadUpdateReconstruct(const fvMesh::readUpdateState stat);

        //- Synchronise non-conformal structure after read-update-reconstruct
        void unconformReadUpdateReconstruct(const fvMesh::readUpdateState stat);

        //- Return the distribution as an FV field for writing
        inline const labelList& cellProc() const
        {
            return cellProc_;
        }

        //- Access the labels of points for each processor
        inline const labelListList& procPointAddressing() const
        {
            validateProcs();
            return procPointAddressing_;
        }

        //- Access the labels of faces for each processor (see notes above)
        inline const labelListList& procFaceAddressing() const
        {
            validateProcs();
            return procFaceAddressing_;
        }

        //- Access the labels of cells for each processor
        inline const labelListList& procCellAddressing() const
        {
            validateProcs();
            return procCellAddressing_;
        }

        //- Access the labels of faces for each processor (see notes above)
        const PtrList<surfaceLabelField::Boundary>&
            procFaceAddressingBf() const;

        //- Write the decomposed meshes and associated data
        void writeComplete(const bool doSets) const;

        //- Write the decomposed meshes and associated data
        void writeProcs(const bool doSets) const;
};


template<>
inline const word& regionName<domainDecomposition>
(
    const domainDecomposition& region
)
{
    return regionName(region.regionName());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
