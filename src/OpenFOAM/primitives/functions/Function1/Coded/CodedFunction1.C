/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2026 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::wordList Foam::Function1s::Coded<Type>::codeKeys
(
    {"code", "codeInclude"}
);

template<class Type>
const Foam::wordList Foam::Function1s::Coded<Type>::codeDictVars
(
    {word::null, word::null}
);

template<class Type>
const Foam::word Foam::Function1s::Coded<Type>::codeOptions
(
    "codedFunction1Options"
);

template<class Type>
const Foam::wordList Foam::Function1s::Coded<Type>::compileFiles
{
    "codedFunction1Template.C"
};

template<class Type>
const Foam::wordList Foam::Function1s::Coded<Type>::copyFiles
{
    "codedFunction1Template.H"
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Coded<Type>::Coded
(
    const word& name,
    const Function1s::unitConversions& units,
    const dictionary& dict
)
:
    Function1<Type>(name),
    codedBase
    (
        dict,
        codeKeys,
        codeDictVars,
        codeOptions,
        compileFiles,
        copyFiles
    ),
    units_(units)
{
    // Set TemplateType filter variables
    setFilterVariable("TemplateType", pTraits<Type>::typeName);

    // Make verbose if debugging
    setFilterVariable("verbose", Foam::name(bool(debug)));

    this->updateLibrary(dict);

    dictionary redirectDict(dict);
    redirectDict.set(codeName(), codeName());

    redirectFunction1Ptr_ =
        Function1<Type>::New(codeName(), unitAny, unitAny, redirectDict);
}


template<class Type>
Foam::Function1s::Coded<Type>::Coded(const Coded<Type>& cf1)
:
    Function1<Type>(cf1),
    codedBase(cf1),
    units_(cf1.units_),
    redirectFunction1Ptr_(cf1.redirectFunction1Ptr_, false)
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
    return
        units_.value.toStandard
        (
            redirectFunction1Ptr_->value
            (
                units_.x.toUser(x)
            )
        );
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
void Foam::Function1s::Coded<Type>::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    codedBase::write(os);
}


// ************************************************************************* //
