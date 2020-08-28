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
    Foam::fvFieldDecomposer

Description
    Finite Volume volume and surface field decomposer.

SourceFiles
    fvFieldDecomposer.C
    fvFieldDecomposerDecomposeFields.C

\*---------------------------------------------------------------------------*/

#ifndef fvFieldDecomposer_H
#define fvFieldDecomposer_H

#include "fvMesh.H"
#include "directFvPatchFieldMapper.H"
#include "surfaceFields.H"

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

        //- Patch field decomposer class
        class patchFieldDecomposer
        :
            public labelList,
            public directFvPatchFieldMapper
        {
            // Private Member Functions

                labelList alignAddressing
                (
                    const labelUList& addressingSlice,
                    const label addressingOffset
                ) const;


        public:

            // Constructors

                //- Construct given addressing
                patchFieldDecomposer
                (
                    const labelUList& addressingSlice,
                    const label addressingOffset
                );
        };


        //- Processor patch field decomposer class. Maps either owner or
        //  neighbour data (no interpolate anymore - processorFvPatchField
        //  holds neighbour data)
        class processorVolPatchFieldDecomposer
        :
            public labelList,
            public directFvPatchFieldMapper
        {
            // Private Member Functions

                labelList alignAddressing
                (
                    const fvMesh& mesh,
                    const labelUList& addressingSlice
                ) const;


        public:

            //- Construct given addressing
            processorVolPatchFieldDecomposer
            (
                const fvMesh& mesh,
                const labelUList& addressingSlice
            );
        };


private:

    // Private Data

        //- Reference to complete mesh
        const fvMesh& completeMesh_;

        //- Reference to processor mesh
        const fvMesh& procMesh_;

        //- Reference to face addressing
        const labelList& faceAddressing_;

        //- Reference to cell addressing
        const labelList& cellAddressing_;

        //- Reference to boundary addressing
        const labelList& boundaryAddressing_;

        //- List of patch field decomposers
        List<patchFieldDecomposer*> patchFieldDecomposerPtrs_;

        List<processorVolPatchFieldDecomposer*>
            processorVolPatchFieldDecomposerPtrs_;


    // Private Member Functions

        //- Helper: map & optionally flip a (face) field
        template<class Type>
        static tmp<Field<Type>> mapField
        (
            const Field<Type>& field,
            const labelUList& mapAndSign,
            const bool applyFlip
        );


public:

    // Constructors

        //- Construct from components
        fvFieldDecomposer
        (
            const fvMesh& completeMesh,
            const fvMesh& procMesh,
            const labelList& faceAddressing,
            const labelList& cellAddressing,
            const labelList& boundaryAddressing
        );

        //- Disallow default bitwise copy construction
        fvFieldDecomposer(const fvFieldDecomposer&) = delete;


    //- Destructor
    ~fvFieldDecomposer();


    // Member Functions

        //- Decompose volume field
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        decomposeField
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const bool allowUnknownPatchFields = false
        ) const;

        //- Decompose surface field
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        decomposeField
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& field
        ) const;

        template<class GeoField>
        void decomposeFields(const PtrList<GeoField>& fields) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvFieldDecomposer&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvFieldDecomposerDecomposeFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
