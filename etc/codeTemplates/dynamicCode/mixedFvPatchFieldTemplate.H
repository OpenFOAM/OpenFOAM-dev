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
    mixed fvPatchField.

    - without state

SourceFiles
    mixedFvPatchFieldTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef mixedFvPatchFieldTemplate_H
#define mixedFvPatchFieldTemplate_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     A templated MixedValueFvPatchField
\*---------------------------------------------------------------------------*/

class ${typeName}MixedValueFvPatch${FieldType}
:
    public mixedFvPatchField<${TemplateType}>
{
public:

    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from patch and internal field
        ${typeName}MixedValueFvPatch${FieldType}
        (
            const fvPatch&,
            const DimensionedField<${TemplateType}, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ${typeName}MixedValueFvPatch${FieldType}
        (
            const fvPatch&,
            const DimensionedField<${TemplateType}, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping a copy onto a new patch
        ${typeName}MixedValueFvPatch${FieldType}
        (
            const ${typeName}MixedValueFvPatch${FieldType}&,
            const fvPatch&,
            const DimensionedField<${TemplateType}, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        ${typeName}MixedValueFvPatch${FieldType}
        (
            const ${typeName}MixedValueFvPatch${FieldType}&
        );

        //- Construct and return a clone
        virtual tmp<fvPatch${FieldType} > clone() const
        {
            return tmp<fvPatch${FieldType} >
            (
                new ${typeName}MixedValueFvPatch${FieldType}(*this)
            );
        }

        //- Copy constructor setting internal field reference
        ${typeName}MixedValueFvPatch${FieldType}
        (
            const ${typeName}MixedValueFvPatch${FieldType}&,
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
                new ${typeName}MixedValueFvPatch${FieldType}(*this, iF)
            );
        }


    //- Destructor
    virtual ~${typeName}MixedValueFvPatch${FieldType}();


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
