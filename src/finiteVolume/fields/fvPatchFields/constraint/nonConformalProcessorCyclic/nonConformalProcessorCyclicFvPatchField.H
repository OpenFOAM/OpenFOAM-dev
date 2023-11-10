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
    Foam::nonConformalProcessorCyclicFvPatchField

Description
    This boundary condition enables processor communication across cyclic
    patches.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nonConformalProcessorCyclic;
    }
    \endverbatim

See also
    Foam::processorCyclicFvPatchField

SourceFiles
    nonConformalProcessorCyclicFvPatchField.C
    nonConformalProcessorCyclicFvPatchFields.H
    nonConformalProcessorCyclicFvPatchFields.C
    nonConformalProcessorCyclicFvPatchFieldsFwd.H

\*---------------------------------------------------------------------------*/

#ifndef nonConformalProcessorCyclicFvPatchField_H
#define nonConformalProcessorCyclicFvPatchField_H

#include "nonConformalProcessorCyclicFvPatch.H"
#include "processorCyclicFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class nonConformalProcessorCyclicFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalProcessorCyclicFvPatchField
:
    public processorCyclicFvPatchField<Type>
{
    // Private Data

        //- Local reference cast into the processor patch
        const nonConformalProcessorCyclicFvPatch& procPatch_;


public:

    //- Runtime type information
    TypeName(nonConformalProcessorCyclicFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        nonConformalProcessorCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        nonConformalProcessorCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        nonConformalProcessorCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nonConformalProcessorCyclicFvPatchField
        //  onto a new patch
        nonConformalProcessorCyclicFvPatchField
        (
            const nonConformalProcessorCyclicFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        nonConformalProcessorCyclicFvPatchField
        (
            const nonConformalProcessorCyclicFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        nonConformalProcessorCyclicFvPatchField
        (
            const nonConformalProcessorCyclicFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new nonConformalProcessorCyclicFvPatchField<Type>(*this, iF)
            );
        }

        //- Construct and return a clone onto a new patch
        //  setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const fvPatch& patch,
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new nonConformalProcessorCyclicFvPatchField<Type>
                (
                    patch,
                    iF,
                    *this
                )
            );
        }


    //- Destructor
    virtual ~nonConformalProcessorCyclicFvPatchField();


    // Member Functions

        // Access

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const
            {
                return procPatch_.transform();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return pTraits<Type>::rank;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonConformalProcessorCyclicFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
