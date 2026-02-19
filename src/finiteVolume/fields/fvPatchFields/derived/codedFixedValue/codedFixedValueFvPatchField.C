/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "codedFixedValueFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::wordList Foam::codedFixedValueFvPatchField<Type>::codeKeys
(
    {"code", "codeInclude", "localCode"}
);

template<class Type>
const Foam::wordList Foam::codedFixedValueFvPatchField<Type>::codeDictVars
(
    {word::null, word::null, word::null}
);

template<class Type>
const Foam::word Foam::codedFixedValueFvPatchField<Type>::codeOptions
(
    "codedFixedValueFvPatchFieldOptions"
);

template<class Type>
const Foam::wordList Foam::codedFixedValueFvPatchField<Type>::compileFiles
{
    "codedFixedValueFvPatchFieldTemplate.C"
};

template<class Type>
const Foam::wordList Foam::codedFixedValueFvPatchField<Type>::copyFiles
{
    "codedFixedValueFvPatchFieldTemplate.H"
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    codedBase
    (
        dict,
        codeKeys,
        codeDictVars,
        codeOptions,
        compileFiles,
        copyFiles
    )
{
    // Set TemplateType and FieldType filter variables
    // (for fvPatchField)
    word fieldType(pTraits<Type>::typeName);

    // Template type for fvPatchField
    setFilterVariable("TemplateType", fieldType);

    // Name for fvPatchField - eg, ScalarField, VectorField, ...
    fieldType[0] = toupper(fieldType[0]);
    setFilterVariable("FieldType", fieldType + "Field");

    // Make verbose if debugging
    setFilterVariable("verbose", Foam::name(bool(debug)));

    // Compile the library containing user-defined fvPatchField
    updateLibrary(dict);
}


template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const codedFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    codedBase(ptf)
{}


template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const codedFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    codedBase(ptf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::fvPatchField<Type>&
Foam::codedFixedValueFvPatchField<Type>::redirectPatchField() const
{
    if (!redirectPatchFieldPtr_.valid())
    {
        OStringStream os;
        writeEntry(os, "type", codeName());

        if (this->overridesConstraint())
        {
            writeEntry(os, "patchType", this->patch().type());
        }

        writeEntry(os, "value", *this);

        IStringStream is(os.str());
        dictionary dict(is);

        redirectPatchFieldPtr_.set
        (
            fvPatchField<Type>::New
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
void Foam::codedFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvPatchField<Type>& fvp = redirectPatchField();

    const_cast<fvPatchField<Type>&>(fvp).updateCoeffs();

    // Copy through value
    this->operator==(fvp);

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    const fvPatchField<Type>& fvp = redirectPatchField();

    const_cast<fvPatchField<Type>&>(fvp).evaluate(commsType);

    fixedValueFvPatchField<Type>::evaluate(commsType);
}


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fixedValueFvPatchField<Type>::write(os);
    codedBase::write(os);
}


// ************************************************************************* //
