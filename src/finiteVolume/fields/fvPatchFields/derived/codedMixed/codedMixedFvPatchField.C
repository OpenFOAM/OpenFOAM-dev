/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
#include "stringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::word Foam::codedMixedFvPatchField<Type>::codeTemplateC
    = "mixedFvPatchFieldTemplate.C";

template<class Type>
const Foam::word Foam::codedMixedFvPatchField<Type>::codeTemplateH
    = "mixedFvPatchFieldTemplate.H";


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
void Foam::codedMixedFvPatchField<Type>::setFieldTemplates
(
    dynamicCode& dynCode
)
{
    word fieldType(pTraits<Type>::typeName);

    // template type for fvPatchField
    dynCode.setFilterVariable("TemplateType", fieldType);

    // Name for fvPatchField - eg, ScalarField, VectorField, ...
    fieldType[0] = toupper(fieldType[0]);
    dynCode.setFilterVariable("FieldType", fieldType + "Field");
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
const Foam::IOdictionary& Foam::codedMixedFvPatchField<Type>::dict() const
{
    const objectRegistry& obr = this->db();

    if (obr.foundObject<IOdictionary>("codeDict"))
    {
        return obr.lookupObject<IOdictionary>("codeDict");
    }
    else
    {
        return obr.store
        (
            new IOdictionary
            (
                IOobject
                (
                    "codeDict",
                    this->db().time().system(),
                    this->db(),
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}


template<class Type>
Foam::dlLibraryTable& Foam::codedMixedFvPatchField<Type>::libs() const
{
    return const_cast<dlLibraryTable&>(this->db().time().libs());
}


template<class Type>
void Foam::codedMixedFvPatchField<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    // take no chances - typeName must be identical to redirectType_
    dynCode.setFilterVariable("typeName", redirectType_);

    // set TemplateType and FieldType filter variables
    // (for fvPatchField)
    setFieldTemplates(dynCode);

    // compile filtered C template
    dynCode.addCompileFile(codeTemplateC);

    // copy filtered H template
    dynCode.addCopyFile(codeTemplateH);


    // debugging: make BC verbose
    //  dynCode.setFilterVariable("verbose", "true");
    //  Info<<"compile " << redirectType_ << " sha1: "
    //      << context.sha1() << endl;

    // define Make/options
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
const Foam::dictionary& Foam::codedMixedFvPatchField<Type>::codeDict()
const
{
    // use system/codeDict or in-line
    return
    (
        dict_.found("code")
      ? dict_
      : this->dict().subDict(redirectType_)
    );
}


template<class Type>
Foam::string Foam::codedMixedFvPatchField<Type>::description() const
{
    return
        "patch "
      + this->patch().name()
      + " on field "
      + this->dimensionedInternalField().name();
}


template<class Type>
void Foam::codedMixedFvPatchField<Type>::clearRedirect() const
{
    // remove instantiation of fvPatchField provided by library
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
    codedBase(),
    dict_(ptf.dict_),
    redirectType_(ptf.redirectType_),
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
    codedBase(),
    dict_(dict),
    redirectType_(dict.lookup("redirectType")),
    redirectPatchFieldPtr_()
{
    updateLibrary(redirectType_);
}


template<class Type>
Foam::codedMixedFvPatchField<Type>::codedMixedFvPatchField
(
    const codedMixedFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    codedBase(),
    dict_(ptf.dict_),
    redirectType_(ptf.redirectType_),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedMixedFvPatchField<Type>::codedMixedFvPatchField
(
    const codedMixedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    codedBase(),
    dict_(ptf.dict_),
    redirectType_(ptf.redirectType_),
    redirectPatchFieldPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::mixedFvPatchField<Type>&
Foam::codedMixedFvPatchField<Type>::redirectPatchField() const
{
    if (!redirectPatchFieldPtr_.valid())
    {
        // Construct a patch
        // Make sure to construct the patchfield with up-to-date value

        // Write the data from the mixed b.c.
        OStringStream os;
        mixedFvPatchField<Type>::write(os);
        IStringStream is(os.str());
        // Construct dictionary from it.
        dictionary dict(is);

        // Override the type to enforce the fvPatchField::New constructor
        // to choose our type
        dict.set("type", redirectType_);

        redirectPatchFieldPtr_.set
        (
            dynamic_cast<mixedFvPatchField<Type>*>
            (
                fvPatchField<Type>::New
                (
                    this->patch(),
                    this->dimensionedInternalField(),
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
    updateLibrary(redirectType_);

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
    updateLibrary(redirectType_);

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
    os.writeKeyword("redirectType") << redirectType_
        << token::END_STATEMENT << nl;

    if (dict_.found("codeInclude"))
    {
        os.writeKeyword("codeInclude")
            << token::HASH << token::BEGIN_BLOCK;

        os.writeQuoted(string(dict_["codeInclude"]), false)
            << token::HASH << token::END_BLOCK
            << token::END_STATEMENT << nl;
    }

    if (dict_.found("localCode"))
    {
        os.writeKeyword("localCode")
            << token::HASH << token::BEGIN_BLOCK;

        os.writeQuoted(string(dict_["localCode"]), false)
            << token::HASH << token::END_BLOCK
            << token::END_STATEMENT << nl;
    }

    if (dict_.found("code"))
    {
        os.writeKeyword("code")
            << token::HASH << token::BEGIN_BLOCK;

        os.writeQuoted(string(dict_["code"]), false)
            << token::HASH << token::END_BLOCK
            << token::END_STATEMENT << nl;
    }

    if (dict_.found("codeOptions"))
    {
        os.writeKeyword("codeOptions")
            << token::HASH << token::BEGIN_BLOCK;

        os.writeQuoted(string(dict_["codeOptions"]), false)
            << token::HASH << token::END_BLOCK
            << token::END_STATEMENT << nl;
    }

    if (dict_.found("codeLibs"))
    {
        os.writeKeyword("codeLibs")
            << token::HASH << token::BEGIN_BLOCK;

        os.writeQuoted(string(dict_["codeLibs"]), false)
            << token::HASH << token::END_BLOCK
            << token::END_STATEMENT << nl;
    }
}


// ************************************************************************* //
