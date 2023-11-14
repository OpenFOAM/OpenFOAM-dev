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
    Foam::fvFieldDecomposer

Description
    Finite Volume volume and surface field decomposer.

SourceFiles
    fvFieldDecomposer.C
    fvFieldDecomposerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fvFieldDecomposer_H
#define fvFieldDecomposer_H

#include "fvMesh.H"
#include "IOobjectList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "forwardFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobjectList;

/*---------------------------------------------------------------------------*\
                      Class fvFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class fvFieldDecomposer
{
public:

    // Public Classes

        //- Patch field decomposer class
        class patchFieldDecomposer
        :
            public labelList,
            public forwardFieldMapper
        {
        public:

            // Constructors

                //- Construct given addressing
                patchFieldDecomposer(const labelUList& addressing);
        };


private:

    // Private Data

        //- Reference to complete mesh
        const fvMesh& completeMesh_;

        //- List of processor meshes
        const PtrList<fvMesh>& procMeshes_;

        //- Reference to face addressing
        const labelListList& faceProcAddressing_;

        //- Reference to cell addressing
        const labelListList& cellProcAddressing_;

        //- Reference to face addressing boundary field
        const PtrList<surfaceLabelField::Boundary>& faceProcAddressingBf_;

        //- List of patch field decomposers
        PtrList<PtrList<patchFieldDecomposer>> patchFieldDecomposers_;


    // Private Member Functions

        //- Convert a processor patch to the corresponding complete patch index
        label completePatchID(const label proci, const label procPatchi) const;

        //- Map cell values to faces
        template<class Type>
        static tmp<Field<Type>> mapCellToFace
        (
            const labelUList& owner,
            const labelUList& neighbour,
            const Field<Type>& field,
            const labelUList& addressing
        );

        //- Map face values to faces
        template<class Type>
        static tmp<Field<Type>> mapFaceToFace
        (
            const Field<Type>& field,
            const labelUList& addressing,
            const bool isFlux
        );

        //- Decompose a volume internal field
        template<class Type>
        PtrList<typename VolField<Type>::Internal>
        decomposeVolInternalField(const IOobject& fieldIoObject) const;

        //- Decompose a volume field
        template<class Type>
        PtrList<VolField<Type>>
        decomposeVolField(const IOobject& fieldIoObject) const;

        //- Decompose a surface field
        template<class Type>
        PtrList<SurfaceField<Type>>
        decomposeFvSurfaceField(const IOobject& fieldIoObject) const;


public:

    // Constructors

        //- Construct from components
        fvFieldDecomposer
        (
            const fvMesh& completeMesh,
            const PtrList<fvMesh>& procMeshes,
            const labelListList& faceProcAddressing,
            const labelListList& cellProcAddressing,
            const PtrList<surfaceLabelField::Boundary>& faceProcAddressingBf
        );

        //- Disallow default bitwise copy construction
        fvFieldDecomposer(const fvFieldDecomposer&) = delete;


    //- Destructor
    ~fvFieldDecomposer();


    // Member Functions

        //- Return whether anything in the object list gets decomposed
        static bool decomposes(const IOobjectList& objects);

        //- Read, decompose and write all volume internal fields
        template<class Type>
        void decomposeVolInternalFields(const IOobjectList& objects);

        //- Read, decompose and write all volume fields
        template<class Type>
        void decomposeVolFields(const IOobjectList& objects);

        //- Read, decompose and write all surface fields
        template<class Type>
        void decomposeFvSurfaceFields(const IOobjectList& objects);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvFieldDecomposer&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvFieldDecomposerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
