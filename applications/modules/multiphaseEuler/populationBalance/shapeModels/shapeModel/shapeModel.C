/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2025 OpenFOAM Foundation
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

#include "shapeModel.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
    defineTypeNameAndDebug(shapeModel, 0);
    defineRunTimeSelectionTable(shapeModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::shapeModel::shapeModel
(
    const populationBalanceModel& popBal
)
:
    popBal_(popBal)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::populationBalance::shapeModel>
Foam::populationBalance::shapeModel::New
(
    const dictionary& dict,
    const populationBalanceModel& popBal
)
{
    const bool haveModelDict = dict.isDict(typeName);

    word modelType;
    const dictionary* modelDictPtr = nullptr;
    if (haveModelDict)
    {
        modelDictPtr = &dict.subDict(typeName);
        modelType = modelDictPtr->lookup<word>("type");
    }
    else
    {
        modelType = dict.lookup<word>(typeName);
        modelDictPtr = &dict.optionalSubDict(modelType + "Coeffs");
    }
    const dictionary& modelDict = *modelDictPtr;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName << " type "
            << modelType << endl << endl
            << "Valid " << typeName << " types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(modelDict, popBal);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalance::shapeModel::~shapeModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::populationBalanceModel&
Foam::populationBalance::shapeModel::popBal() const
{
    return popBal_;
}


void Foam::populationBalance::shapeModel::solve()
{}


void Foam::populationBalance::shapeModel::correct()
{}


void Foam::populationBalance::shapeModel::addCoalescence
(
    const volScalarField::Internal& Su,
    const label i,
    const label j,
    const label k
)
{}


void Foam::populationBalance::shapeModel::addBreakup
(
    const volScalarField::Internal& Su,
    const label i,
    const label j
)
{}


void Foam::populationBalance::shapeModel::reset()
{}


// ************************************************************************* //
