/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::mappedValueAndPatchInternalValueFvPatchField

Description
    This boundary condition maps the boundary and internal values from a
    neighbouring patch field to the boundary and internal values of this patch.

Usage
    \table
        Property     | Description                | Required | Default value
        field        | name of field to be mapped | no       | this field name
        setAverage   | set the average value?     | no       | yes if average \\
                                                               is specified, \\
                                                               no otherwise
        average      | average value to apply     | if setAverage is true |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            mappedValueAndPatchInternalValue;
        field           T;
        average         300;
        value           uniform 300;
    }
    \endverbatim

    This boundary condition will usually be applied to a patch which is of
    mappedPatchBase type, and which holds all the necessary mapping
    information. It can also create its own mapping data which overrides that
    in the mapped patch, or so that it can be applied to a non-mapped patch.
    This is triggered by the presence of controls relating to mappedPatchBase
    (i.e., neighbourRegion, neighbourPatch, etc ...).

See also
    Foam::mappedPatchBase
    Foam::mappedPolyPatch
    Foam::mappedFvPatch
    Foam::mappedValueFvPatchField

SourceFiles
    mappedValueAndPatchInternalValue.C

\*---------------------------------------------------------------------------*/

#ifndef mappedValueAndPatchInternalValueFvPatchField_H
#define mappedValueAndPatchInternalValueFvPatchField_H

#include "mappedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class mappedValueAndPatchInternalValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mappedValueAndPatchInternalValueFvPatchField
:
    public mappedValueFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("mappedValueAndPatchInternalValue");


    // Constructors

        //- Construct from patch and internal field
        mappedValueAndPatchInternalValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mappedValueAndPatchInternalValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mappedValueAndPatchInternalValue onto a new patch
        mappedValueAndPatchInternalValueFvPatchField
        (
            const mappedValueAndPatchInternalValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        mappedValueAndPatchInternalValueFvPatchField
        (
            const mappedValueAndPatchInternalValueFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        mappedValueAndPatchInternalValueFvPatchField
        (
            const mappedValueAndPatchInternalValueFvPatchField<Type>&,
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
                new mappedValueAndPatchInternalValueFvPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedValueAndPatchInternalValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
