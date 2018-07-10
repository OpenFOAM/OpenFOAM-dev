/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    Foam::jumpCyclicAMIFvPatchField

Description
    This boundary condition provides a base class that enforces a cyclic
    condition with a specified 'jump' (or offset) between a pair of boundaries,
    whereby communication between the patches is performed using an arbitrary
    mesh interface (AMI) interpolation.

See also
    Foam::cyclicAMIFvPatchField

SourceFiles
    jumpCyclicAMIFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef jumpCyclicAMIFvPatchField_H
#define jumpCyclicAMIFvPatchField_H

#include "cyclicAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class jumpCyclicAMIFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class jumpCyclicAMIFvPatchField
:
    public cyclicAMIFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("jumpCyclicAMI");


    // Constructors

        //- Construct from patch and internal field
        jumpCyclicAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        jumpCyclicAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given jumpCyclicAMIFvPatchField onto a
        //  new patch
        jumpCyclicAMIFvPatchField
        (
            const jumpCyclicAMIFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        jumpCyclicAMIFvPatchField
        (
            const jumpCyclicAMIFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        jumpCyclicAMIFvPatchField
        (
            const jumpCyclicAMIFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );


    // Member functions

        // Access

            //- Return the interface type
            virtual const word& interfaceFieldType() const
            {
                return cyclicAMIFvPatchField<Type>::type();
            }

            //- Return the "jump" across the patch as a "half" field
            virtual tmp<Field<Type>> jump() const = 0;


        // Evaluation functions

            //- Return neighbour coupled given internal cell data
            tmp<Field<Type>> patchNeighbourField() const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                scalarField& result,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>&,
                const Field<Type>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const;
};


//- Update result field based on interface functionality
template<>
void jumpCyclicAMIFvPatchField<scalar>::updateInterfaceMatrix
(
    scalarField& result,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "jumpCyclicAMIFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
