/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

Description
    Template for use with dynamic code generation of a
    fixedValue fvPatchField.

    - without state

SourceFiles
    fixedValueFvPatchFieldTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValueFvPatchFieldTemplate_H
#define fixedValueFvPatchFieldTemplate_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         A templated FixedValueFvPatch
\*---------------------------------------------------------------------------*/

class ${typeName}FixedValueFvPatch${FieldType}
:
    public fixedValueFvPatchField<${TemplateType}>
{
public:

    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from patch and internal field
        ${typeName}FixedValueFvPatch${FieldType}
        (
            const fvPatch&,
            const DimensionedField<${TemplateType}, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ${typeName}FixedValueFvPatch${FieldType}
        (
            const fvPatch&,
            const DimensionedField<${TemplateType}, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping a copy onto a new patch
        ${typeName}FixedValueFvPatch${FieldType}
        (
            const ${typeName}FixedValueFvPatch${FieldType}&,
            const fvPatch&,
            const DimensionedField<${TemplateType}, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        ${typeName}FixedValueFvPatch${FieldType}
        (
            const ${typeName}FixedValueFvPatch${FieldType}&
        ) = delete;

        //- Copy constructor setting internal field reference
        ${typeName}FixedValueFvPatch${FieldType}
        (
            const ${typeName}FixedValueFvPatch${FieldType}&,
            const DimensionedField<${TemplateType}, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatch${FieldType} > clone
        (
            const DimensionedField<${TemplateType}, volMesh>& iF
        ) const
        {
            return tmp<fvPatch${FieldType} >
            (
                new ${typeName}FixedValueFvPatch${FieldType}(*this, iF)
            );
        }


    //- Destructor
    virtual ~${typeName}FixedValueFvPatch${FieldType}();


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
