/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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

#include "codedFixedValuePointPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "pointPatchFieldMapper.H"
#include "pointFields.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::codedFixedValuePointPatchField<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    // Take no chances - typeName must be identical to codeName()
    dynCode.setFilterVariable("typeName", codeName());

    // Set TemplateType and FieldType filter variables
    // (for pointPatchField)
    word fieldType(pTraits<Type>::typeName);

    // Template type for pointPatchField
    dynCode.setFilterVariable("TemplateType", fieldType);

    // Name for pointPatchField - eg, ScalarField, VectorField, ...
    fieldType[0] = toupper(fieldType[0]);
    dynCode.setFilterVariable("FieldType", fieldType + "Field");

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC);

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH);

    // Debugging: make BC verbose
    if (debug)
    {
        // Debugging: make BC verbose
        dynCode.setFilterVariable("verbose", "true");
        Info<<"compile " << codeName() << " sha1: "
            << context.sha1() << endl;
    }

    // Define Make/options
    dynCode.setMakeOptions
    (
        "EXE_INC = -g \\\n"
        "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
      + context.options()
      + "\n\nLIB_LIBS = \\\n"
      + "    -lOpenFOAM \\\n"
      + "    -lfiniteVolume \\\n"
      + context.libs()
    );
}


template<class Type>
void Foam::codedFixedValuePointPatchField<Type>::clearRedirect() const
{
    // Remove instantiation of pointPatchField provided by library
    redirectPatchFieldPtr_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::codedFixedValuePointPatchField<Type>::codedFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(p, iF),
    CodedBase<codedFixedValuePointPatchFieldBase>(),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedFixedValuePointPatchField<Type>::codedFixedValuePointPatchField
(
    const codedFixedValuePointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<Type>(ptf, p, iF, mapper),
    CodedBase<codedFixedValuePointPatchFieldBase>(ptf),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedFixedValuePointPatchField<Type>::codedFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<Type>(p, iF, dict),
    CodedBase<codedFixedValuePointPatchFieldBase>(dict),
    redirectPatchFieldPtr_()
{
    updateLibrary();
}


template<class Type>
Foam::codedFixedValuePointPatchField<Type>::codedFixedValuePointPatchField
(
    const codedFixedValuePointPatchField<Type>& ptf
)
:
    fixedValuePointPatchField<Type>(ptf),
    CodedBase<codedFixedValuePointPatchFieldBase>(ptf),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedFixedValuePointPatchField<Type>::codedFixedValuePointPatchField
(
    const codedFixedValuePointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(ptf, iF),
    CodedBase<codedFixedValuePointPatchFieldBase>(ptf),
    redirectPatchFieldPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::pointPatchField<Type>&
Foam::codedFixedValuePointPatchField<Type>::redirectPatchField() const
{
    if (!redirectPatchFieldPtr_.valid())
    {
        OStringStream os;
        writeEntry(os, "type", codeName());
        writeEntry(os, "value", static_cast<const Field<Type>&>(*this));
        IStringStream is(os.str());
        dictionary dict(is);

        redirectPatchFieldPtr_.set
        (
            pointPatchField<Type>::New
            (
                this->patch(),
                this->internalField(),
                dict
            ).ptr()
        );
    }

    return redirectPatchFieldPtr_();
}


template<class Type>
void Foam::codedFixedValuePointPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Make sure library containing user-defined pointPatchField is up-to-date
    updateLibrary();

    const pointPatchField<Type>& fvp = redirectPatchField();

    const_cast<pointPatchField<Type>&>(fvp).updateCoeffs();

    // Copy through value
    this->operator==(fvp);

    fixedValuePointPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::codedFixedValuePointPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Make sure library containing user-defined pointPatchField is up-to-date
    updateLibrary();

    const pointPatchField<Type>& fvp = redirectPatchField();

    const_cast<pointPatchField<Type>&>(fvp).evaluate(commsType);

    fixedValuePointPatchField<Type>::evaluate(commsType);
}


template<class Type>
void Foam::codedFixedValuePointPatchField<Type>::write(Ostream& os) const
{
    fixedValuePointPatchField<Type>::write(os);
    writeCode(os);
}


// ************************************************************************* //
