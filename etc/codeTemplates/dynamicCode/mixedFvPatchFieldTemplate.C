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

\*---------------------------------------------------------------------------*/

#include "mixedFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode
${localCode}
//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ${SHA1sum}
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void ${typeName}_${SHA1sum}(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatch${FieldType},
    ${typeName}MixedValueFvPatch${FieldType}
);


const char* const ${typeName}MixedValueFvPatch${FieldType}::SHA1sum =
    "${SHA1sum}";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

${typeName}MixedValueFvPatch${FieldType}::
${typeName}MixedValueFvPatch${FieldType}
(
    const fvPatch& p,
    const DimensionedField<${TemplateType}, volMesh>& iF
)
:
    mixedFvPatchField<${TemplateType}>(p, iF)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " from patch/DimensionedField\n";
    }
}


${typeName}MixedValueFvPatch${FieldType}::
${typeName}MixedValueFvPatch${FieldType}
(
    const ${typeName}MixedValueFvPatch${FieldType}& ptf,
    const fvPatch& p,
    const DimensionedField<${TemplateType}, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<${TemplateType}>(ptf, p, iF, mapper)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " from patch/DimensionedField/mapper\n";
    }
}


${typeName}MixedValueFvPatch${FieldType}::
${typeName}MixedValueFvPatch${FieldType}
(
    const fvPatch& p,
    const DimensionedField<${TemplateType}, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<${TemplateType}>(p, iF, dict)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " from patch/dictionary\n";
    }
}


${typeName}MixedValueFvPatch${FieldType}::
${typeName}MixedValueFvPatch${FieldType}
(
    const ${typeName}MixedValueFvPatch${FieldType}& ptf,
    const DimensionedField<${TemplateType}, volMesh>& iF
)
:
    mixedFvPatchField<${TemplateType}>(ptf, iF)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum} "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

${typeName}MixedValueFvPatch${FieldType}::
~${typeName}MixedValueFvPatch${FieldType}()
{
    if (${verbose:-false})
    {
        Info<<"destroy ${typeName} sha1: ${SHA1sum}\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ${typeName}MixedValueFvPatch${FieldType}::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (${verbose:-false})
    {
        Info<<"updateCoeffs ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin code
    ${code}
//}}} end code

    this->mixedFvPatchField<${TemplateType}>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
