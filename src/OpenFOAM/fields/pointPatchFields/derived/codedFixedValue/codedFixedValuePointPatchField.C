/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
#include "stringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::word Foam::codedFixedValuePointPatchField<Type>::codeTemplateC
    = "fixedValuePointPatchFieldTemplate.C";

template<class Type>
const Foam::word Foam::codedFixedValuePointPatchField<Type>::codeTemplateH
    = "fixedValuePointPatchFieldTemplate.H";


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
void Foam::codedFixedValuePointPatchField<Type>::setFieldTemplates
(
    dynamicCode& dynCode
)
{
    word fieldType(pTraits<Type>::typeName);

    // Template type for pointPatchField
    dynCode.setFilterVariable("TemplateType", fieldType);

    // Name for pointPatchField - eg, ScalarField, VectorField, ...
    fieldType[0] = toupper(fieldType[0]);
    dynCode.setFilterVariable("FieldType", fieldType + "Field");
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
const Foam::IOdictionary& Foam::codedFixedValuePointPatchField<Type>::dict()
const
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
Foam::dlLibraryTable& Foam::codedFixedValuePointPatchField<Type>::libs() const
{
    return const_cast<dlLibraryTable&>(this->db().time().libs());
}


template<class Type>
void Foam::codedFixedValuePointPatchField<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    // Take no chances - typeName must be identical to name_
    dynCode.setFilterVariable("typeName", name_);

    // Set TemplateType and FieldType filter variables
    // (for pointPatchField)
    setFieldTemplates(dynCode);

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC);

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH);


    // Debugging: make BC verbose
    //   dynCode.setFilterVariable("verbose", "true");
    //   Info<<"compile " << name_ << " sha1: "
    //       << context.sha1() << endl;

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
const Foam::dictionary& Foam::codedFixedValuePointPatchField<Type>::codeDict()
const
{
    // Use system/codeDict or in-line
    return
    (
        dict_.found("code")
      ? dict_
      : this->dict().subDict(name_)
    );
}


template<class Type>
Foam::string Foam::codedFixedValuePointPatchField<Type>::description() const
{
    return
        "patch "
      + this->patch().name()
      + " on field "
      + this->internalField().name();
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
    codedBase(),
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
    codedBase(),
    dict_(ptf.dict_),
    name_(ptf.name_),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedFixedValuePointPatchField<Type>::codedFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedValuePointPatchField<Type>(p, iF, dict, valueRequired),
    codedBase(),
    dict_(dict),
    name_
    (
        dict.found("redirectType")
      ? dict.lookup("redirectType")
      : dict.lookup("name")
    ),
    redirectPatchFieldPtr_()
{
    updateLibrary(name_);
}


template<class Type>
Foam::codedFixedValuePointPatchField<Type>::codedFixedValuePointPatchField
(
    const codedFixedValuePointPatchField<Type>& ptf
)
:
    fixedValuePointPatchField<Type>(ptf),
    codedBase(),
    dict_(ptf.dict_),
    name_(ptf.name_),
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
    codedBase(),
    dict_(ptf.dict_),
    name_(ptf.name_),
    redirectPatchFieldPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::pointPatchField<Type>&
Foam::codedFixedValuePointPatchField<Type>::redirectPatchField() const
{
    if (!redirectPatchFieldPtr_.valid())
    {
        // Construct a patch
        // Make sure to construct the patchfield with up-to-date value

        OStringStream os;
        os.writeKeyword("type") << name_ << token::END_STATEMENT
            << nl;
        static_cast<const Field<Type>&>(*this).writeEntry("value", os);
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
    updateLibrary(name_);

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
    updateLibrary(name_);

    const pointPatchField<Type>& fvp = redirectPatchField();

    const_cast<pointPatchField<Type>&>(fvp).evaluate(commsType);

    fixedValuePointPatchField<Type>::evaluate(commsType);
}


template<class Type>
void Foam::codedFixedValuePointPatchField<Type>::write(Ostream& os) const
{
    fixedValuePointPatchField<Type>::write(os);
    os.writeKeyword("name") << name_
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
