/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::fixedValueInletOutletFvPatchField

Description
    This boundary condition sets a fixed value. When the flow direction is
    inwards this acts exactly like a fixed value condition. In the presence of
    outflow, however, this condition approximates the fixed value constraint in
    advective terms by fixing the gradient instead.

    This condition is not likely to be used on its own. It is more suitable as
    a base class for conditions that need to specify the value of a field even
    when the flow reverses.

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | Name of the flux field  | no          | phi
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedValueInletOutlet;
        phi             phi;
        value           0;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    fixedValueInletOutletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValueInletOutletFvPatchField_H
#define fixedValueInletOutletFvPatchField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fixedValueInletOutletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedValueInletOutletFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private Data

        //- Name of the flux field
        const word phiName_;


public:

    //- Runtime type information
    TypeName("fixedValueInletOutlet");


    // Constructors

        //- Construct from patch and internal field
        fixedValueInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedValueInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping given fixedValueInletOutletFvPatchField
        //  onto a new patch
        fixedValueInletOutletFvPatchField
        (
            const fixedValueInletOutletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );

        //- Disallow copy without setting internal field reference
        fixedValueInletOutletFvPatchField
        (
            const fixedValueInletOutletFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        fixedValueInletOutletFvPatchField
        (
            const fixedValueInletOutletFvPatchField<Type>&,
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
                new fixedValueInletOutletFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fixedValueInletOutletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
