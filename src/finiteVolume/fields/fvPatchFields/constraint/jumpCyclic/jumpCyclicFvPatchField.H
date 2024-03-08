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
    Foam::jumpCyclicFvPatchField

Description
    This boundary condition provides a base class for cyclic conditions with a
    specified "jump" (or offset) between the sides

See also
    Foam::cyclicFvPatchField

SourceFiles
    jumpCyclicFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef jumpCyclicFvPatchField_H
#define jumpCyclicFvPatchField_H

#include "cyclicFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class jumpCyclicFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class jumpCyclicFvPatchField
:
    public cyclicFvPatchField<Type>
{
    // Private Member Functions

        //- Update result field based on interface functionality
        template<class ... Cmpt>
        void updateInterfaceMatrixCmpts
        (
            Field<Type>& result,
            const Field<Type>& psi,
            const scalarField& coeffs,
            const Pstream::commsTypes comms,
            const Cmpt ... cmpt
        ) const;


public:

    //- Runtime type information
    TypeName("jumpCyclic");


    // Constructors

        //- Construct from patch and internal field
        jumpCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        jumpCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given jumpCyclicFvPatchField onto a new patch
        jumpCyclicFvPatchField
        (
            const jumpCyclicFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        jumpCyclicFvPatchField(const jumpCyclicFvPatchField<Type>&) = delete;

        //- Copy constructor setting internal field reference
        jumpCyclicFvPatchField
        (
            const jumpCyclicFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );


    // Member Functions

        // Access

            //- Return the interface type
            virtual const word& interfaceFieldType() const
            {
                return cyclicFvPatchField<Type>::type();
            }

            //- Return the "jump" in value applied to the opposite field when
            //  transforming to this side
            virtual tmp<Field<Type>> jump() const = 0;


        // Evaluation functions

            //- Return neighbour coupled given internal cell data
            virtual tmp<Field<Type>> patchNeighbourField
            (
                const Pstream::commsTypes
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                scalarField& result,
                const scalarField& psi,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes comms
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>& result,
                const Field<Type>& psi,
                const scalarField& coeffs,
                const Pstream::commsTypes comms
            ) const;
};


//- Update result field based on interface functionality
template<>
void jumpCyclicFvPatchField<scalar>::updateInterfaceMatrix
(
    scalarField& result,
    const scalarField& psi,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes comms
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "jumpCyclicFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
