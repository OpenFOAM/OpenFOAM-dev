/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
#include "labelList.H"
#include "SLList.H"
#include "PtrList.H"
#include "point.H"
#include "Time.H"
#include "volFields.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class domainDecomposition Declaration
\*---------------------------------------------------------------------------*/

class domainDecomposition
:
    public fvMesh
{
    // Private Data

        //- Optional: points at the facesInstance
        autoPtr<pointIOField> facesInstancePointsPtr_;

        //- Number of processors in decomposition
        label nProcs_;

        //- Is the decomposition data to be distributed for each processor
        bool distributed_;

        //- Processor label for each cell
        labelList cellToProc_;

        //- Labels of points for each processor
        labelListList procPointAddressing_;

        //- Labels of faces for each processor
        // Note: Face turning index is stored as the sign on addressing
        // Only the processor boundary faces are affected: if the sign of the
        // index is negative, the processor face is the reverse of the
        // original face. In order to do this properly, all face
        // indices will be incremented by 1 and the decremented as
        // necessary to avoid the problem of face number zero having no
        // sign.
        List<DynamicList<label>> procFaceAddressing_;

        //- Labels of cells for each processor
        labelListList procCellAddressing_;

        //- Sizes for processor mesh patches
        // Excludes inter-processor boundaries
        labelListList procPatchSize_;

        //- Start indices for processor patches
        // Excludes inter-processor boundaries
        labelListList procPatchStartIndex_;


        // Per inter-processor patch information

        //- Neighbour processor ID for inter-processor boundaries
        labelListList procNeighbourProcessors_;

        //- Sizes for inter-processor patches
        labelListList procProcessorPatchSize_;

        //- Start indices (in procFaceAddressing_) for inter-processor patches
        labelListList procProcessorPatchStartIndex_;

        //- Sub patch IDs for inter-processor patches
        List<labelListList> procProcessorPatchSubPatchIDs_;

        //- Sub patch sizes for inter-processor patches
        List<labelListList> procProcessorPatchSubPatchStarts_;

    // Private Member Functions

        void distributeCells(const fileName& dictFile);

        //- Mark all elements with value or -2 if occur twice
        static void mark
        (
            const labelList& zoneElems,
            const label zoneI,
            labelList& elementToZone
        );

        //- Append single element to list
        static void append(labelList&, const label);

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
        void processInterCyclics
        (
            const polyBoundaryMesh& patches,
            List<DynamicList<DynamicList<label>>>& interPatchFaces,
            List<Map<label>>& procNbrToInterPatch,
            List<labelListList>& subPatchIDs,
            List<labelListList>& subPatchStarts,
            bool owner,
            BinaryOp bop
        ) const;


public:

    // Constructors

        //- Construct from IOobject and decomposition dictionary name
        domainDecomposition
        (
            const IOobject& io,
            const fileName& dictFile
        );


    //- Destructor
    ~domainDecomposition();


    // Member Functions

        //- Number of processor in decomposition
        label nProcs() const
        {
            return nProcs_;
        }

        //- Is the decomposition data to be distributed for each processor
        bool distributed() const
        {
            return distributed_;
        }

        //- Decompose mesh.
        void decomposeMesh(const fileName& dict);

        //- Write decomposition
        bool writeDecomposition(const bool decomposeSets);

        //- Cell-processor decomposition labels
        const labelList& cellToProc() const
        {
            return cellToProc_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "domainDecompositionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
