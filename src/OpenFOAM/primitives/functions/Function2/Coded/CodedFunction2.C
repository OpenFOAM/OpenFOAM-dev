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

#include "CodedFunction2.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::Function2s::Coded<Type>::codeKeys() const
{
    return
    {
        "code",
        "codeInclude"
    };
}


template<class Type>
void Foam::Function2s::Coded<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    dynCode.setFilterVariable("typeName", codeName());

    // Set TemplateType filter variables
    dynCode.setFilterVariable("TemplateType", pTraits<Type>::typeName);

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC("codedFunction2"));

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH("codedFunction2"));

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
void Foam::Function2s::Coded<Type>::clearRedirect() const
{
    // Remove instantiation of Function2 provided by library
    redirectFunction2Ptr_.clear();
}


template<class Type>
Foam::autoPtr<Foam::Function2<Type>>
Foam::Function2s::Coded<Type>::compileNew()
{
    this->updateLibrary();

    dictionary redirectDict(codeDict());
    redirectDict.set(codeName(), codeName());

    return Function2<Type>::New(codeName(), redirectDict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2s::Coded<Type>::Coded
(
    const word& name,
    const dictionary& dict
)
:
    Function2<Type>(name),
    codedBase(name, dict)
{
    redirectFunction2Ptr_ = compileNew();
}



template<class Type>
Foam::Function2s::Coded<Type>::Coded(const Coded<Type>& cf1)
:
    Function2<Type>(cf1),
    codedBase(cf1)
{
    redirectFunction2Ptr_ = compileNew();
}


template<class Type>
Foam::tmp<Foam::Function2<Type>> Foam::Function2s::Coded<Type>::clone() const
{
    return tmp<Function2<Type>>(new Coded<Type>(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2s::Coded<Type>::~Coded()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function2s::Coded<Type>::value
(
    const scalarField& x,
    const scalarField& y
) const
{
    return redirectFunction2Ptr_->value(x, y);
}


template<class Type>
void Foam::Function2s::Coded<Type>::write(Ostream& os) const
{
    writeCode(os);
}


// ************************************************************************* //
