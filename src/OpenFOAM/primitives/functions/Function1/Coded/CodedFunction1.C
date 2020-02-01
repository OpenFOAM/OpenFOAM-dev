/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "CodedFunction1.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"

// * * * * * * * * * * * * Private Static Data Members * * * * * * * * * * * //

template<class Type>
const Foam::wordList Foam::Function1s::Coded<Type>::codeKeys_ =
{
    "code",
    "codeInclude"
};


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::word Foam::Function1s::Coded<Type>::codeTemplateC =
    "Function1Template.C";

template<class Type>
const Foam::word Foam::Function1s::Coded<Type>::codeTemplateH =
    "Function1Template.H";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::Function1s::Coded<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    dynCode.setFilterVariable("typeName", name_);

    // Set TemplateType filter variables
    dynCode.setFilterVariable("TemplateType", pTraits<Type>::typeName);

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC);

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH);

    // debugging: make verbose
    //  dynCode.setFilterVariable("verbose", "true");
    //  Info<<"compile " << name_ << " sha1: "
    //      << context.sha1() << endl;

    // define Make/options
    dynCode.setMakeOptions
        (
            "EXE_INC = -g \\\n"
            + context.options()
            + "\n\nLIB_LIBS = \\\n"
            + "    -lOpenFOAM \\\n"
            + context.libs()
        );
}


template<class Type>
Foam::string Foam::Function1s::Coded<Type>::description() const
{
    return Function1<Type>::typeName_() + (" " + name_);
}


template<class Type>
void Foam::Function1s::Coded<Type>::clearRedirect() const
{
    // Remove instantiation of Function1 provided by library
    redirectFunction1Ptr_.clear();
}


template<class Type>
const Foam::dictionary& Foam::Function1s::Coded<Type>::codeDict()
const
{
    return dict_;
}


template<class Type>
const Foam::wordList& Foam::Function1s::Coded<Type>::codeKeys() const
{
    return codeKeys_;
}


template<class Type>
Foam::autoPtr<Foam::Function1<Type>>
Foam::Function1s::Coded<Type>::compileAndLink()
{
    updateLibrary(name_);

    dictionary redirectDict(dict_);
    redirectDict.set(name_, name_);

    return Function1<Type>::New(name_, redirectDict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Coded<Type>::Coded
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<Type>(entryName),
    dict_(dict),
    name_
    (
        dict.found("redirectType")
      ? dict.lookup("redirectType")
      : dict.lookup("name")
    ),
    redirectFunction1Ptr_(compileAndLink())
{}



template<class Type>
Foam::Function1s::Coded<Type>::Coded(const Coded<Type>& cf1)
:
    Function1<Type>(cf1),
    codedBase(),
    dict_(cf1.dict_),
    name_(cf1.name_),
    redirectFunction1Ptr_(compileAndLink())
{}


template<class Type>
Foam::tmp<Foam::Function1<Type>> Foam::Function1s::Coded<Type>::clone() const
{
    return tmp<Function1<Type>>(new Coded<Type>(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Coded<Type>::~Coded()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1s::Coded<Type>::value
(
    const scalarField& x
) const
{
    return redirectFunction1Ptr_->value(x);
}


template<class Type>
inline Type Foam::Function1s::Coded<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
    return pTraits<Type>::zero;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1s::Coded<Type>::integrate
(
    const scalarField& x1,
    const scalarField& x2
) const
{
    NotImplemented;
    return tmp<Field<Type>>();
}


template<class Type>
void Foam::Function1s::Coded<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os  << token::END_STATEMENT << nl;
    writeEntry(os, "name", name_);

    if (dict_.found("codeInclude"))
    {
        writeKeyword(os, "codeInclude");
        os.write(verbatimString(dict_["codeInclude"]))
            << token::END_STATEMENT << nl;
    }

    if (dict_.found("code"))
    {
        writeKeyword(os, "code");
        os.write(verbatimString(dict_["code"]))
            << token::END_STATEMENT << nl;
    }
}


// ************************************************************************* //
