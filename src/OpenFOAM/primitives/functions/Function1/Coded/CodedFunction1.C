/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::Function1s::Coded<Type>::codeKeys() const
{
    return
    {
        "code",
        "codeInclude"
    };
}


template<class Type>
void Foam::Function1s::Coded<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    dynCode.setFilterVariable("typeName", codeName());

    // Set TemplateType filter variables
    dynCode.setFilterVariable("TemplateType", pTraits<Type>::typeName);

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC("codedFunction1"));

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH("codedFunction1"));

    // Debugging: make verbose
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
      + context.options()
      + "\n\nLIB_LIBS = \\\n"
      + "    -lOpenFOAM \\\n"
      + context.libs()
    );
}


template<class Type>
void Foam::Function1s::Coded<Type>::clearRedirect() const
{
    // Remove instantiation of Function1 provided by library
    redirectFunction1Ptr_.clear();
}


template<class Type>
Foam::autoPtr<Foam::Function1<Type>>
Foam::Function1s::Coded<Type>::compileNew()
{
    this->updateLibrary();

    dictionary redirectDict(codeDict());
    redirectDict.set(codeName(), codeName());

    return Function1<Type>::New(codeName(), redirectDict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Coded<Type>::Coded
(
    const word& name,
    const dictionary& dict
)
:
    Function1<Type>(name),
    codedBase(dict)
{
    redirectFunction1Ptr_ = compileNew();
}



template<class Type>
Foam::Function1s::Coded<Type>::Coded(const Coded<Type>& cf1)
:
    Function1<Type>(cf1),
    codedBase(cf1)
{
    redirectFunction1Ptr_ = compileNew();
}


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
inline Type Foam::Function1s::Coded<Type>::integral
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
    return pTraits<Type>::zero;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1s::Coded<Type>::integral
(
    const scalarField& x1,
    const scalarField& x2
) const
{
    NotImplemented;
    return tmp<Field<Type>>();
}


template<class Type>
void Foam::Function1s::Coded<Type>::write(Ostream& os) const
{
    writeCode(os);
}


// ************************************************************************* //
