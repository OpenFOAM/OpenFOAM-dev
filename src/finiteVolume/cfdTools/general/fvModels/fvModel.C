/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "fvModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvModel, 0);
    defineRunTimeSelectionTable(fvModel, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fvModel::addSupType
(
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{}


template<class Type>
void Foam::fvModel::addSupType
(
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{}


template<class Type>
void Foam::fvModel::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvModel::fvModel
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    modelType_(modelType),
    mesh_(mesh)
{
    Info<< incrIndent << indent << "Name: " << name_
        << endl << decrIndent;
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fvModel> Foam::fvModel::New
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
{
    const word modelType(dict.lookup("type"));

    Info<< indent
        << "Selecting finite volume model type " << modelType << endl;

    if
    (
        !dictionaryConstructorTablePtr_
     || dictionaryConstructorTablePtr_->find(modelType)
        == dictionaryConstructorTablePtr_->end()
    )
    {
        if
        (
           !libs.open
            (
                dict,
                "libs",
                dictionaryConstructorTablePtr_
            )
        )
        {
            libs.open("lib" + modelType.remove(':') + ".so", false);
        }

        if (!dictionaryConstructorTablePtr_)
        {
            FatalErrorInFunction
                << "Unknown model type "
                << modelType << nl << nl
                << "Table of fvModels is empty"
                << exit(FatalError);
        }
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown fvModel " << modelType << nl << nl
            << "Valid fvModels are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<fvModel>
    (
        cstrIter()(name, modelType, mesh, dict)
    );
}


Foam::fvModel::~fvModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fvModel::addSupFields() const
{
    return wordList::null();
}


bool Foam::fvModel::addsSupToField(const word& fieldName) const
{
    return findIndex(addSupFields(), fieldName) != -1;
}


Foam::scalar Foam::fvModel::maxDeltaT() const
{
    return vGreat;
}


void Foam::fvModel::addSup(fvMatrix<scalar>& eqn) const
{}


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_FIELD_SUP, fvModel)


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_RHO_FIELD_SUP, fvModel)


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP, fvModel)


void Foam::fvModel::preUpdateMesh()
{}


void Foam::fvModel::correct()
{}


bool Foam::fvModel::read(const dictionary& dict)
{
    return true;
}


bool Foam::fvModel::write(const bool write) const
{
    return true;
}


// ************************************************************************* //
