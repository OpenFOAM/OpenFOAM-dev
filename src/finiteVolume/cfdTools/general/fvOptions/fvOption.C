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

#include "fvOption.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(option, 0);
        defineRunTimeSelectionTable(option, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fv::option::addSupType
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{}


template<class Type>
void Foam::fv::option::addSupType
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{}


template<class Type>
void Foam::fv::option::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{}


template<class Type>
void Foam::fv::option::constrainType
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{}


template<class Type>
void Foam::fv::option::correctType(VolField<Type>& field) const
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::option::option
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    name_(name),
    modelType_(modelType),
    mesh_(mesh),
    dict_(dict),
    coeffs_(dict.optionalSubDict(modelType + "Coeffs"))
{
    Info<< incrIndent << indent << "Source: " << name_ << endl << decrIndent;
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fv::option> Foam::fv::option::New
(
    const word& name,
    const dictionary& coeffs,
    const fvMesh& mesh
)
{
    word modelType(coeffs.lookup("type"));

    Info<< indent
        << "Selecting finite volume options model type " << modelType << endl;

    libs.open
    (
        coeffs,
        "libs",
        dictionaryConstructorTablePtr_
    );

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown Model type " << modelType << nl << nl
            << "Valid model types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<option>(cstrIter()(name, modelType, coeffs, mesh));
}


Foam::fv::option::~option()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::option::addSupFields() const
{
    return wordList::null();
}


Foam::wordList Foam::fv::option::constrainedFields() const
{
    return wordList::null();
}


Foam::wordList Foam::fv::option::correctedFields() const
{
    return wordList::null();
}


bool Foam::fv::option::addsSupToField(const word& fieldName) const
{
    return findIndex(addSupFields(), fieldName) != -1;
}


bool Foam::fv::option::constrainsField(const word& fieldName) const
{
    return findIndex(constrainedFields(), fieldName) != -1;
}


bool Foam::fv::option::correctsField(const word& fieldName) const
{
    return findIndex(correctedFields(), fieldName) != -1;
}


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_OPTION_ADD_SUP, option);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_OPTION_ADD_RHO_SUP, option);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_OPTION_ADD_ALPHA_RHO_SUP, option);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_OPTION_CONSTRAIN, option);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_OPTION_CORRECT, option);


void Foam::fv::option::updateMesh(const mapPolyMesh& mpm)
{}


bool Foam::fv::option::movePoints()
{
    return true;
}


// ************************************************************************* //
