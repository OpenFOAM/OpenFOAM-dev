/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

/*---------------------------------------------------------------------------*\
                     Class domainDecomposition Declaration
\*---------------------------------------------------------------------------*/

class domainDecomposition
{
    // Private Data

        //- Run times
        const processorRunTimes& runTimes_;

        //- Region name
        const word regionName_;

        //- The complete mesh
        autoPtr<fvMesh> completeMesh_;

        //- The processor meshes
        PtrList<fvMesh> procMeshes_;


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

            //- ...
            List<List<DynamicList<label>>>
                nonConformalProcFaceAddressingBf() const;

            // Map the polyFacesBf from the processor to the complete meshes and
            // unconform the complete mesh
            void unconformComplete();

            // Map the polyFacesBf from the complete to the processor meshes and
            // unconform the processor meshes
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
        void readComplete(const bool stitch = true);

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
            const word& regionName
        );


    //- Destructor
    virtual ~domainDecomposition();


    // Member Functions

        //- Access the run times
        inline const processorRunTimes& procRunTimes() const
        {
            return runTimes_;
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
        bool readDecompose(const bool doSets);

        //- Read in the processor meshes. Read the complete mesh if it is
        //  available and up to date relative to the processor meshes, or
        //  reconstruct if not. Return whether or not reconstruction happened.
        bool readReconstruct(const bool doSets);

        //- Read-update for decomposition
        fvMesh::readUpdateState readUpdateDecompose();

        //- Read-update for reconstruction
        fvMesh::readUpdateState readUpdateReconstruct();

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
