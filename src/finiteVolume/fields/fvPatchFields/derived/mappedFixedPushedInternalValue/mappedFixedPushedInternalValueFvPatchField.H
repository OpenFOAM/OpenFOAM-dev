/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::mappedFixedPushedInternalValueFvPatchField

Description
    This boundary condition maps the boundary values of a neighbour patch
    field to the boundary and internal cell values of *this.

Usage
    \table
        Property     | Description             | Required    | Default value
        fieldName    | name of field to be mapped | no       | this field name
        setAverage   | flag to activate setting of average value | yes |
        average      | average value to apply if \c setAverage = yes | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            mappedFixedPushedInternalValue;
        fieldName       T;
        setAverage      no;
        average         0;
        value           uniform 0;
    }
    \endverbatim

Note
    This boundary condition can only be applied to patches that are of
    the \c mappedPolyPatch type.

See also
    Foam::mappedPatchBase
    Foam::mappedPolyPatch
    Foam::mappedFvPatch
    Foam::mappedFixedValueFvPatchField

SourceFiles
    mappedFixedPushedInternalValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedFixedPushedInternalValueFvPatchField_H
#define mappedFixedPushedInternalValueFvPatchField_H

#include "mappedFixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class mappedFixedPushedInternalValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mappedFixedPushedInternalValueFvPatchField
:
    public mappedFixedValueFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("mappedFixedPushedInternalValue");


    // Constructors

        //- Construct from patch and internal field
        mappedFixedPushedInternalValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mappedFixedPushedInternalValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        //  mappedFixedPushedInternalValueFvPatchField  onto a new patch
        mappedFixedPushedInternalValueFvPatchField
        (
            const mappedFixedPushedInternalValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        mappedFixedPushedInternalValueFvPatchField
        (
            const mappedFixedPushedInternalValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new mappedFixedPushedInternalValueFvPatchField<Type>
                (
                    *this
                )
            );
        }

        //- Copy constructor setting internal field reference
        mappedFixedPushedInternalValueFvPatchField
        (
            const mappedFixedPushedInternalValueFvPatchField<Type>&,
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
                new mappedFixedPushedInternalValueFvPatchField<Type>
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
    #include "mappedFixedPushedInternalValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
