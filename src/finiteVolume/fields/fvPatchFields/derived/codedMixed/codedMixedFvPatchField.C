/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "codedMixedFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::codedMixedFvPatchField<Type>::codeKeys() const
{
    return
    {
        "code",
        "codeInclude",
        "localCode"
    };
}


template<class Type>
void Foam::codedMixedFvPatchField<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    dynCode.setFilterVariable("typeName", codeName());

    // Set TemplateType and FieldType filter variables
    // (for fvPatchField)
    word fieldType(pTraits<Type>::typeName);

    // Template type for fvPatchField
    dynCode.setFilterVariable("TemplateType", fieldType);

    // Name for fvPatchField - eg, ScalarField, VectorField, ...
    fieldType[0] = toupper(fieldType[0]);
    dynCode.setFilterVariable("FieldType", fieldType + "Field");

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC("codedMixedFvPatchField"));

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH("codedMixedFvPatchField"));

    // Debugging: make BC verbose
    if (debug)
    {
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
void Foam::codedMixedFvPatchField<Type>::clearRedirect() const
{
    // Remove instantiation of fvPatchField provided by library
    redirectPatchFieldPtr_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::codedMixedFvPatchField<Type>::codedMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    codedBase(),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedMixedFvPatchField<Type>::codedMixedFvPatchField
(
    const codedMixedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    codedBase(ptf),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedMixedFvPatchField<Type>::codedMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF, dict),
    codedBase(dict),
    redirectPatchFieldPtr_()
{
    updateLibrary();
}


template<class Type>
Foam::codedMixedFvPatchField<Type>::codedMixedFvPatchField
(
    const codedMixedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    codedBase(ptf),
    redirectPatchFieldPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::mixedFvPatchField<Type>&
Foam::codedMixedFvPatchField<Type>::redirectPatchField() const
{
    if (!redirectPatchFieldPtr_.valid())
    {
        OStringStream os;
        mixedFvPatchField<Type>::write(os);
        IStringStream is(os.str());
        dictionary dict(is);

        // Override the type to enforce the fvPatchField::New constructor
        // to choose our type
        dict.set("type", codeName());

        redirectPatchFieldPtr_.set
        (
            dynamic_cast<mixedFvPatchField<Type>*>
            (
                fvPatchField<Type>::New
                (
                    this->patch(),
                    this->internalField(),
                    dict
                ).ptr()
            )
        );
    }

    return redirectPatchFieldPtr_();
}


template<class Type>
void Foam::codedMixedFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Make sure library containing user-defined fvPatchField is up-to-date
    updateLibrary();

    const mixedFvPatchField<Type>& fvp = redirectPatchField();

    const_cast<mixedFvPatchField<Type>&>(fvp).updateCoeffs();

    // Copy through coefficients
    this->refValue() = fvp.refValue();
    this->refGrad() = fvp.refGrad();
    this->valueFraction() = fvp.valueFraction();

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::codedMixedFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Make sure library containing user-defined fvPatchField is up-to-date
    updateLibrary();

    const mixedFvPatchField<Type>& fvp = redirectPatchField();

    // - updates the value of fvp (though not used)
    // - resets the updated() flag
    const_cast<mixedFvPatchField<Type>&>(fvp).evaluate(commsType);

    // Update the value (using the coefficients) locally
    mixedFvPatchField<Type>::evaluate(commsType);
}


template<class Type>
void Foam::codedMixedFvPatchField<Type>::write(Ostream& os) const
{
    mixedFvPatchField<Type>::write(os);
    writeCode(os);
}


// ************************************************************************* //
