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
    fixedValue pointPatchField.

    - without state

SourceFiles
    fixedValuePointPatchFieldTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValuePointPatchFieldTemplate_H
#define fixedValuePointPatchFieldTemplate_H

#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                A templated FixedValuePointPatchField
\*---------------------------------------------------------------------------*/

class ${typeName}FixedValuePointPatch${FieldType}
:
    public fixedValuePointPatchField<${TemplateType}>
{
public:

    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from patch and internal field
        ${typeName}FixedValuePointPatch${FieldType}
        (
            const pointPatch&,
            const DimensionedField<${TemplateType}, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ${typeName}FixedValuePointPatch${FieldType}
        (
            const pointPatch&,
            const DimensionedField<${TemplateType}, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping a copy onto a new patch
        ${typeName}FixedValuePointPatch${FieldType}
        (
            const ${typeName}FixedValuePointPatch${FieldType}&,
            const pointPatch&,
            const DimensionedField<${TemplateType}, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Copy constructor
        ${typeName}FixedValuePointPatch${FieldType}
        (
            const ${typeName}FixedValuePointPatch${FieldType}&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<${TemplateType}>> clone() const
        {
            return autoPtr<pointPatchField<${TemplateType}>>
            (
                new ${typeName}FixedValuePointPatch${FieldType}(*this)
            );
        }

        //- Copy constructor setting internal field reference
        ${typeName}FixedValuePointPatch${FieldType}
        (
            const ${typeName}FixedValuePointPatch${FieldType}&,
            const DimensionedField<${TemplateType}, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<${TemplateType}>> clone
        (
            const DimensionedField<${TemplateType}, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<${TemplateType}>>
            (
                new ${typeName}FixedValuePointPatch${FieldType}(*this, iF)
            );
        }


    //- Destructor
    virtual ~${typeName}FixedValuePointPatch${FieldType}();


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
