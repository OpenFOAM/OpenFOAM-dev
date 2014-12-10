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

#include "codedFixedValueFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::word Foam::codedFixedValueFvPatchField<Type>::codeTemplateC
    = "fixedValueFvPatchFieldTemplate.C";

template<class Type>
const Foam::word Foam::codedFixedValueFvPatchField<Type>::codeTemplateH
    = "fixedValueFvPatchFieldTemplate.H";


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::setFieldTemplates
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
const Foam::IOdictionary& Foam::codedFixedValueFvPatchField<Type>::dict() const
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
Foam::dlLibraryTable& Foam::codedFixedValueFvPatchField<Type>::libs() const
{
    return const_cast<dlLibraryTable&>(this->db().time().libs());
}


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::prepare
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
const Foam::dictionary& Foam::codedFixedValueFvPatchField<Type>::codeDict()
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
Foam::string Foam::codedFixedValueFvPatchField<Type>::description() const
{
    return
        "patch "
      + this->patch().name()
      + " on field "
      + this->dimensionedInternalField().name();
}


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::clearRedirect() const
{
    // remove instantiation of fvPatchField provided by library
    redirectPatchFieldPtr_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    codedBase(),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const codedFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    codedBase(),
    dict_(ptf.dict_),
    redirectType_(ptf.redirectType_),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    codedBase(),
    dict_(dict),
    redirectType_(dict.lookup("redirectType")),
    redirectPatchFieldPtr_()
{
    updateLibrary(redirectType_);
}


template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const codedFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    codedBase(),
    dict_(ptf.dict_),
    redirectType_(ptf.redirectType_),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedFixedValueFvPatchField<Type>::codedFixedValueFvPatchField
(
    const codedFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    codedBase(),
    dict_(ptf.dict_),
    redirectType_(ptf.redirectType_),
    redirectPatchFieldPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::fvPatchField<Type>&
Foam::codedFixedValueFvPatchField<Type>::redirectPatchField() const
{
    if (!redirectPatchFieldPtr_.valid())
    {
        // Construct a patch
        // Make sure to construct the patchfield with up-to-date value

        OStringStream os;
        os.writeKeyword("type") << redirectType_ << token::END_STATEMENT
            << nl;
        static_cast<const Field<Type>&>(*this).writeEntry("value", os);
        IStringStream is(os.str());
        dictionary dict(is);

        redirectPatchFieldPtr_.set
        (
            fvPatchField<Type>::New
            (
                this->patch(),
                this->dimensionedInternalField(),
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

    // Make sure library containing user-defined fvPatchField is up-to-date
    updateLibrary(redirectType_);

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
    // Make sure library containing user-defined fvPatchField is up-to-date
    updateLibrary(redirectType_);

    const fvPatchField<Type>& fvp = redirectPatchField();

    const_cast<fvPatchField<Type>&>(fvp).evaluate(commsType);

    fixedValueFvPatchField<Type>::evaluate(commsType);
}


template<class Type>
void Foam::codedFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fixedValueFvPatchField<Type>::write(os);
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
