/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "Coded_DimensionedFieldFunction.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::wordList Foam::DimensionedFieldFunctions::Coded<Type>::codeKeys
(
    {"evaluate", "update", "codeInclude"}
);

template<class Type>
const Foam::wordList Foam::DimensionedFieldFunctions::Coded<Type>::codeDictVars
(
    {word::null, word::null, word::null}
);

template<class Type>
const Foam::word Foam::DimensionedFieldFunctions::Coded<Type>::codeOptions
(
    "codedDimensionedFieldFunctionOptions"
);

template<class Type>
const Foam::wordList Foam::DimensionedFieldFunctions::Coded<Type>::compileFiles
{
    "codedDimensionedFieldFunctionTemplate.C"
};

template<class Type>
const Foam::wordList Foam::DimensionedFieldFunctions::Coded<Type>::copyFiles
{
    "codedDimensionedFieldFunctionTemplate.H"
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class DimensionedFieldType>
Foam::DimensionedFieldFunctions::Coded<DimensionedFieldType>::
Coded
(
    const dictionary& dict,
    DimensionedFieldType& field
)
:
    DimensionedFieldFunction<DimensionedFieldType>(dict, field),
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
    // Set variable substitutions
    varSubstitutions().set
    (
        {
            {"DimensionedFieldType", DimensionedFieldType::typeName},
            {
                "DimensionedFieldTypeName",
                codedName(DimensionedFieldType::typeName)
            },
            {"verbose", Foam::name(bool(debug))}
        }
    );

    this->updateLibrary(dict);

    dictionary redirectDict(dict);
    redirectDict.set("type", codeName());

    redirectFunctionPtr_ = DimensionedFieldFunction<DimensionedFieldType>::New
    (
        redirectDict,
        field
    );

    evaluate();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class DimensionedFieldType>
void Foam::DimensionedFieldFunctions::Coded<DimensionedFieldType>::
evaluate()
{
    redirectFunctionPtr_->evaluate();
}


template<class DimensionedFieldType>
void Foam::DimensionedFieldFunctions::Coded<DimensionedFieldType>::write
(
    Ostream& os
) const
{
    codedBase::write(os);
}


// ************************************************************************* //
