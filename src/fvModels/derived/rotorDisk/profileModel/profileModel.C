/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "profileModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(profileModel, 0);
    defineRunTimeSelectionTable(profileModel, dictionary);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::profileModel::readFromFile() const
{
    return fName_ != fileName::null;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::profileModel::profileModel(const dictionary& dict, const word& name)
:
    dict_(dict),
    name_(name),
    fName_(fileName::null)
{
    dict.readIfPresent("file", fName_);
}

// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::profileModel::~profileModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::profileModel::name() const
{
    return name_;
}


Foam::autoPtr<Foam::profileModel> Foam::profileModel::New
(
    const dictionary& dict
)
{
    const word& modelName(dict.dictName());

    const word modelType(dict.lookup("type"));

    Info<< "    - creating " << modelType << " profile " << modelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown profile model type " << modelType
            << nl << nl
            << "Valid model types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<profileModel>(cstrIter()(dict, modelName));
}


// ************************************************************************* //
