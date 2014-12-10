/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "distributionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace distributionModels
    {
        defineTypeNameAndDebug(distributionModel, 0);
        defineRunTimeSelectionTable(distributionModel, dictionary);
    }
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::distributionModels::distributionModel::check() const
{
    if (minValue() < 0)
    {
        FatalErrorIn("distributionModels::distributionModel::check() const")
            << type() << "distribution: Minimum value must be greater than "
            << "zero." << nl << "Supplied minValue = " << minValue()
            << abort(FatalError);
    }

    if (maxValue() < minValue())
    {
        FatalErrorIn("distributionModels::distributionModel::check() const")
            << type() << "distribution: Maximum value is smaller than the "
            << "minimum value:" << nl << "    maxValue = " << maxValue()
            << ", minValue = " << minValue()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::distributionModel::distributionModel
(
    const word& name,
    const dictionary& dict,
    cachedRandom& rndGen
)
:
    distributionModelDict_(dict.subDict(name + "Distribution")),
    rndGen_(rndGen)
{}


Foam::distributionModels::distributionModel::distributionModel
(
    const distributionModel& p
)
:
    distributionModelDict_(p.distributionModelDict_),
    rndGen_(p.rndGen_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributionModels::distributionModel::~distributionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::distributionModel::sample() const
{
    notImplemented
    (
        "Foam::scalar "
        "Foam::distributionModels::distributionModel::sample() const"
    );
    return 0.0;
}


Foam::scalar Foam::distributionModels::distributionModel::minValue() const
{
    notImplemented
    (
        "Foam::scalar "
        "Foam::distributionModels::distributionModel::minValue() const"
    );
    return 0.0;
}


Foam::scalar Foam::distributionModels::distributionModel::maxValue() const
{
    notImplemented
    (
        "Foam::scalar "
        "Foam::distributionModels::distributionModel::maxValue() const"
    );
    return 0.0;
}


Foam::scalar Foam::distributionModels::distributionModel::meanValue() const
{
    notImplemented
    (
        "Foam::scalar "
        "Foam::distributionModels::distributionModel::meanValue() const"
    );
    return 0.0;
}


// ************************************************************************* //
