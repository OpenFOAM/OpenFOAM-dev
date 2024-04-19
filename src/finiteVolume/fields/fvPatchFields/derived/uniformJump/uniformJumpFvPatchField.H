/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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
    Foam::uniformJumpFvPatchField

Description
    This boundary condition provides a jump condition, using the \c cyclic
    condition as a base. The jump is specified as the difference between the
    neighbour patch and the owner patch value (i.e., neighbour minus owner).
    The jump value is uniform across the patches, and can vary in time.

Usage
    \table
        Property     | Description             | Required    | Default value
        patchType    | underlying patch type (should be \c cyclic) | yes |
        jumpTable    | jump value              | yes         |
    \endtable

    Example of the boundary condition specification with a constant jump:
    \verbatim
    <patchName>
    {
        type            uniformJump;
        patchType       cyclic;
        jumpTable       constant 10;
    }
    \endverbatim

    Example of the boundary condition specification with a time-varying jump:
    \verbatim
    <patchName>
    {
        type            uniformJump;
        patchType       cyclic;
        jumpTable
        {
            type            table;
            values
            (
                (0 0)
                (0.1 10)
                (0.9 10)
                (1 0)
            );
        }
    }
    \endverbatim

See also
    Foam::fixedJumpFvPatchField
    Foam::Function1s

SourceFiles
    uniformJumpFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformJumpFvPatchField_H
#define uniformJumpFvPatchField_H

#include "fixedJumpFvPatchField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class uniformJumpFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class uniformJumpFvPatchField
:
    public fixedJumpFvPatchField<Type>
{
protected:

    // Protected data

        //- The "jump" table
        autoPtr<Function1<Type>> jumpTable_;


public:

    //- Runtime type information
    TypeName("uniformJump");


    // Constructors

        //- Construct from patch, internal field and dictionary
        uniformJumpFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given uniformJumpFvPatchField onto a
        //  new patch
        uniformJumpFvPatchField
        (
            const uniformJumpFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        uniformJumpFvPatchField(const uniformJumpFvPatchField<Type>&) = delete;

        //- Copy constructor setting internal field reference
        uniformJumpFvPatchField
        (
            const uniformJumpFvPatchField<Type>&,
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
                new uniformJumpFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "uniformJumpFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
