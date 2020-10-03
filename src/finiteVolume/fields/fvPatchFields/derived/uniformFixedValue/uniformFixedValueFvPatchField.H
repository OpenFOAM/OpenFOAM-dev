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
    Foam::uniformFixedValueFvPatchField

Description
    This boundary condition provides a uniform fixed value condition.

Usage
    \table
        Property     | Description             | Required    | Default value
        uniformValue | uniform value           | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            uniformFixedValue;
        uniformValue    constant 0.2;
    }
    \endverbatim

Note
    The uniformValue entry is a Function1 type, able to describe time
    varying functions.  The example above gives the usage for supplying a
    constant value.

See also
    Foam::Function1s
    Foam::fixedValueFvPatchField

SourceFiles
    uniformFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformFixedValueFvPatchField_H
#define uniformFixedValueFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class uniformFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class uniformFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private Data

        autoPtr<Function1<Type>> uniformValue_;


public:

    //- Runtime type information
    TypeName("uniformFixedValue");


    // Constructors

        //- Construct from patch and internal field
        uniformFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        uniformFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>& fld
        );

        //- Construct from patch, internal field and dictionary
        uniformFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given uniformFixedValueFvPatchField
        //  onto a new patch
        uniformFixedValueFvPatchField
        (
            const uniformFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        uniformFixedValueFvPatchField
        (
            const uniformFixedValueFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        uniformFixedValueFvPatchField
        (
            const uniformFixedValueFvPatchField<Type>&,
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
                new uniformFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "uniformFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
