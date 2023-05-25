/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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

#include "sinteringModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace shapeModels
{
    defineTypeNameAndDebug(sinteringModel, 0);
    defineRunTimeSelectionTable(sinteringModel, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::shapeModels::sinteringModel::sinteringModel
(
    const dictionary& dict,
    const fractal& fractalShape
)
:
    fractal_(fractalShape)
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diameterModels::shapeModels::sinteringModel>
Foam::diameterModels::shapeModels::sinteringModel::New
(
    const dictionary& dict,
    const fractal& fractalShape
)
{
    word sinteringModelType(dict.lookup("sinteringModel"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(sinteringModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown sinteringModel type "
            << sinteringModelType << endl << endl
            << "Valid sinteringModel types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict, fractalShape);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::shapeModels::sinteringModel::~sinteringModel()
{}


// ************************************************************************* //
