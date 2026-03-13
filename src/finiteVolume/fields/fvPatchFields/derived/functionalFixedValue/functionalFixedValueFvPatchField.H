/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::functionalFixedValueFvPatchField

Description
    This boundary condition extrapolates field to the patch using the near-cell
    values and adjusts the distribution to match the specified, optionally
    time-varying, mean value.

Usage
    \table
        Property     | Description             | Required    | Default value
        meanValue    | mean value Function1    | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            functionalFixedValue;
        meanValue       1.0;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    functionalFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef functionalFixedValueFvPatchField_H
#define functionalFixedValueFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "UDimensionedFvPatchFields.H"
#include "DimensionedFvPatchFieldsFwd.H"
#include "DimensionedFieldFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class functionalFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class functionalFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private Data

        UDimensionedFvPatchField<Type> dimensionedValue_;

        autoPtr<DimensionedFieldFunction<DimensionedFvPatchField<Type>>>
            funcPtr_;

public:

    //- Runtime type information
    TypeName("functionalFixedValue");


    // Constructors

        //- Construct from patch, internal field and dictionary
        functionalFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given functionalFixedValueFvPatchField
        //  onto a new patch
        functionalFixedValueFvPatchField
        (
            const functionalFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, fvMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        functionalFixedValueFvPatchField
        (
            const functionalFixedValueFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        functionalFixedValueFvPatchField
        (
            const functionalFixedValueFvPatchField<Type>&,
            const DimensionedField<Type, fvMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, fvMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new functionalFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchField<Type>&, const fieldMapper&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchField<Type>&);


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using fixedValueFvPatchField<Type>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "functionalFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
