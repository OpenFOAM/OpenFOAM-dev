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
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class domainDecomposition Declaration
\*---------------------------------------------------------------------------*/

class domainDecomposition
{
    // Private Data

        //- Reference to the complete mesh
        const fvMesh& mesh_;

        //- Number of processors in decomposition
        const label nProcs_;

        //- Optional: points at the facesInstance
        autoPtr<pointIOField> facesInstancePointsPtr_;

        //- Is the decomposition data to be distributed for each processor
        bool distributed_;


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
            List<DynamicList<label>> procFaceAddressing_;

            //- Labels of cells for each processor
            labelListList procCellAddressing_;


        //- Processor times
        PtrList<Time> procRunTimes_;

        //- Processor meshes
        PtrList<fvMesh> procMeshes_;



    // Private Member Functions

        //- Call the decomposition method and return the processor index that
        //  each cell is being distributed to
        labelList distributeCells();

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
            List<Map<label>>&,
            List<DynamicList<DynamicList<label>>>&
        ) const;

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

        //- Validate that the decomposition has been generated or read
        void validate() const;


public:

    //- Runtime type information
    TypeName("domainDecomposition");


    // Constructors

        //- Construct from reference to fvMesh
        domainDecomposition(const fvMesh& mesh);


    //- Destructor
    virtual ~domainDecomposition();


    // Member Functions

        //- Access the global mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return the number of processors in the decomposition
        label nProcs() const
        {
            return nProcs_;
        }

        //- Is the decomposition data to be distributed for each processor?
        bool distributed() const
        {
            return distributed_;
        }

        //- Access the labels of points for each processor
        const labelListList& procPointAddressing() const
        {
            validate();
            return procPointAddressing_;
        }

        //- Access the labels of faces for each processor (see notes above)
        const List<DynamicList<label>>& procFaceAddressing() const
        {
            validate();
            return procFaceAddressing_;
        }

        //- Access the labels of cells for each processor
        const labelListList& procCellAddressing() const
        {
            validate();
            return procCellAddressing_;
        }

        //- Return the processor meshes
        const PtrList<fvMesh>& procMeshes() const
        {
            validate();
            return procMeshes_;
        }

        //- Generate the decomposition
        void decompose();

        //- Reset the time
        void setTime(const instant&, const label newIndex);

        //- Read additional (faceInstance) points, if necessary
        void readPoints();

        //- Read in the decomposition addressing
        void readAddressing();

        //- Read in the meshes and the decomposition
        void read();

        //- Read update the meshes and the decomposition
        fvMesh::readUpdateState readUpdate();

        //- Return the distribution as an FV field for writing
        labelList cellToProc() const;

        //- Read additional (faceInstance) points, if necessary
        void writePoints() const;

        //- Write the decomposition addressing
        void writeAddressing() const;

        //- Write the decomposed meshes and associated data
        void write(const bool decomposeSets) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
