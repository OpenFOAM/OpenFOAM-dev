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
    Foam::fvFieldReconstructor

Description
    Finite volume reconstructor for volume and surface fields.

SourceFiles
    fvFieldReconstructor.C
    fvFieldReconstructorTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fvFieldReconstructor_H
#define fvFieldReconstructor_H

#include "PtrList.H"
#include "surfaceFields.H"
#include "IOobjectList.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fvFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class fvFieldReconstructor
{
    // Private Data

        //- Reconstructed mesh reference
        const fvMesh& completeMesh_;

        //- List of processor meshes
        const PtrList<fvMesh>& procMeshes_;

        //- List of processor face addressing lists
        const labelListList& faceProcAddressing_;

        //- List of processor cell addressing lists
        const labelListList& cellProcAddressing_;

        //- Boundary field of face addressing
        const PtrList<surfaceLabelField::Boundary>& faceProcAddressingBf_;


    // Private Member Functions

        //- Return whether anything in the object list gets reconstructed
        template<class FieldType>
        static bool reconstructs
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        );

        //- Convert a processor patch to the corresponding complete patch index
        label completePatchID(const label proci, const label procPatchi) const;

        //- ...
        template<class Type>
        static void rmapFaceToFace
        (
            Field<Type>& toField,
            const Field<Type>& fromField,
            const labelUList& addressing,
            const bool isFlux
        );

        //- Read and reconstruct a volume internal field
        template<class Type>
        tmp<DimensionedField<Type, volMesh>>
        reconstructVolInternalField(const IOobject& fieldIoObject) const;

        //- Read and reconstruct a volume field
        template<class Type>
        tmp<VolField<Type>>
        reconstructVolField(const IOobject& fieldIoObject) const;

        //- Read and reconstruct a surface field
        template<class Type>
        tmp<SurfaceField<Type>>
        reconstructFvSurfaceField(const IOobject& fieldIoObject) const;


public:

    // Constructors

        //- Construct from components
        fvFieldReconstructor
        (
            const fvMesh& mesh,
            const PtrList<fvMesh>& procMeshes,
            const labelListList& faceProcAddressing,
            const labelListList& cellProcAddressing,
            const PtrList<surfaceLabelField::Boundary>& faceProcAddressingBf
        );

        //- Disallow default bitwise copy construction
        fvFieldReconstructor(const fvFieldReconstructor&) = delete;


    // Member Functions

        //- Return whether anything in the object list gets reconstructed
        static bool reconstructs
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        );

        //- Read, reconstruct and write all/selected volume internal fields
        template<class Type>
        void reconstructVolInternalFields
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        );

        //- Read, reconstruct and write all/selected volume fields
        template<class Type>
        void reconstructVolFields
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        );

        //- Read, reconstruct and write all/selected surface fields
        template<class Type>
        void reconstructFvSurfaceFields
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvFieldReconstructor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvFieldReconstructorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
