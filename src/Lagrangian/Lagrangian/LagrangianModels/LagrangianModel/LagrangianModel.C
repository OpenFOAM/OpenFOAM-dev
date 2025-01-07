/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LagrangianModel, 0);
    defineRunTimeSelectionTable(LagrangianModel, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::LagrangianModel::addSupType
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubSubField<Type>& field,
    LagrangianEqn<Type>& eqn
) const
{}


template<class Type>
void Foam::LagrangianModel::addSupType
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& m,
    const LagrangianSubSubField<Type>& field,
    LagrangianEqn<Type>& eqn
) const
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LagrangianModel::LagrangianModel
(
    const word& name,
    const LagrangianMesh& mesh
)
:
    name_(name),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::LagrangianModel> Foam::LagrangianModel::New
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict
)
{
    const word type(modelDict.lookup("type"));

    Info<< "Selecting " << typeName
        << " with name " << name
        << " of type " << type << endl;

    if
    (
        !dictionaryConstructorTablePtr_
     || dictionaryConstructorTablePtr_->find(type)
     == dictionaryConstructorTablePtr_->end()
    )
    {
        if (!libs.open(modelDict, "libs", dictionaryConstructorTablePtr_))
        {
            libs.open("lib" + type.remove(':') + ".so", false);
        }

        if (!dictionaryConstructorTablePtr_)
        {
            FatalErrorInFunction
                << "Unknown " << typeName << " type "
                << type << nl << nl
                << "Table of " << typeName << "s is empty"
                << exit(FatalError);
        }
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(modelDict)
            << "Unknown " << typeName << " " << type << nl << nl
            << "Valid " << typeName << "s are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return
        cstrIter()
        (
            name,
            mesh,
            modelDict.optionalSubDict(type + "Coeffs"),
            stateDict(name, mesh)
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LagrangianModel::~LagrangianModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::wordList Foam::LagrangianModel::addSupFields() const
{
    return wordList::null();
}


bool Foam::LagrangianModel::addsSupToField(const word& fieldName) const
{
    return findIndex(addSupFields(), fieldName) != -1;
}


void Foam::LagrangianModel::postConstruct()
{}


void Foam::LagrangianModel::correct()
{}


void Foam::LagrangianModel::preModify
(
    const LagrangianMesh& mesh,
    DynamicList<elementModification>& elementModifications
) const
{}


Foam::LagrangianSubMesh Foam::LagrangianModel::modify
(
    LagrangianMesh& mesh,
    const LagrangianSubMesh& modifiedMesh
) const
{
    return mesh_.subNone();
}


void Foam::LagrangianModel::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{}


void Foam::LagrangianModel::addSup
(
    const LagrangianSubScalarField& deltaT,
    LagrangianSubScalarField& S
) const
{}


FOR_ALL_FIELD_TYPES(IMPLEMENT_LAGRANGIAN_MODEL_ADD_FIELD_SUP, LagrangianModel)


FOR_ALL_FIELD_TYPES(IMPLEMENT_LAGRANGIAN_MODEL_ADD_M_FIELD_SUP, LagrangianModel)


void Foam::LagrangianModel::topoChange(const polyTopoChangeMap&)
{}


void Foam::LagrangianModel::mapMesh(const polyMeshMap&)
{}


void Foam::LagrangianModel::distribute(const polyDistributionMap&)
{}


bool Foam::LagrangianModel::read(const dictionary& dict)
{
    return true;
}


bool Foam::LagrangianModel::write(const bool write) const
{
    return write;
}


// ************************************************************************* //
