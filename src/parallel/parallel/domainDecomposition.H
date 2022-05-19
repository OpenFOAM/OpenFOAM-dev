/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    decomposeMesh.C

\*---------------------------------------------------------------------------*/

#ifndef domainDecomposition_H
#define domainDecomposition_H

#include "fvMesh.H"
#include "processorRunTimes.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

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


        // Processor mesh to complete mesh addressing

            //- Labels of points for each processor
            labelListList procPointAddressing_;

            //- Labels of faces for each processor
            // Note: Face turning index is stored as the sign on addressing
            // Only the processor boundary faces are affected: if the sign of
            // the index is negative, the processor face is the reverse of the
            // original face. In order to do this properly, all face
            // indices will be incremented by 1 and the decremented as
            // necessary to avoid the problem of face number zero having no
            // sign.
            labelListList procFaceAddressing_;

            //- Labels of cells for each processor
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
            template<class BinaryOp>
            inline void processInterCyclics
            (
                const labelList& cellToProc,
                const polyBoundaryMesh& patches,
                List<DynamicList<DynamicList<label>>>& interPatchFaces,
                List<Map<label>>& procNbrToInterPatch,
                List<labelListList>& subPatchIDs,
                List<labelListList>& subPatchStarts,
                bool owner,
                BinaryOp bop
            ) const;


        //- ...
        void decomposePoints();

        //- ...
        void reconstructPoints();

        //- ...
        bool completeConformal() const;

        //- ...
        bool procsConformal() const;

        //- If the complete mesh is non-conformal and the processor meshes are
        //  not, then map the non-conformal addressing across and unconform the
        //  processor meshes. And vice versa.
        void unconform();

        //- ...
        label compareInstances(const fileName& a, const fileName& b) const;

        //- Validate that the complete mesh has been generated or read
        void validateComplete() const;

        //- Validate that the processor meshes have been generated or read
        void validateProcs() const;


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

        //- Access the global mesh
        const fvMesh& completeMesh() const
        {
            validateComplete();
            return completeMesh_();
        }

        //- Access the processor meshes
        const PtrList<fvMesh>& procMeshes() const
        {
            validateProcs();
            return procMeshes_;
        }

        //- Return the number of processors in the decomposition
        label nProcs() const
        {
            return runTimes_.nProcs();
        }

        //- ...
        void readComplete();

        //- ...
        void readProcs();

        //- ...
        void readAddressing();

        //- ...
        fvMesh::readUpdateState readUpdate();

        //- Decompose the complete mesh to create the processor meshes and
        //  populate the addressing
        void decompose();

        //- Reconstruct the processor meshes to create the complete mesh and
        //  populate the addressing
        void reconstruct();

        //- Return the distribution as an FV field for writing
        labelList cellToProc() const;

        //- Access the labels of points for each processor
        const labelListList& procPointAddressing() const
        {
            validateProcs();
            return procPointAddressing_;
        }

        //- Access the labels of faces for each processor (see notes above)
        const labelListList& procFaceAddressing() const
        {
            validateProcs();
            return procFaceAddressing_;
        }

        //- Access the labels of cells for each processor
        const labelListList& procCellAddressing() const
        {
            validateProcs();
            return procCellAddressing_;
        }

        //- Access the labels of faces for each processor (see notes above)
        const PtrList<surfaceLabelField::Boundary>&
            procFaceAddressingBf() const;

        //- Write the decomposition addressing
        void writeAddressing() const;

        //- Write the decomposed meshes and associated data
        void writeComplete(const bool doSets) const;

        //- Write the decomposed meshes and associated data
        void writeProcs(const bool doSets) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
