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

#include "ejectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ejectionModel, 0);
defineRunTimeSelectionTable(ejectionModel, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void ejectionModel::addToEjectedMass(const scalar dMass)
{
    ejectedMass_ += dMass;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ejectionModel::ejectionModel(surfaceFilmRegionModel& film)
:
    filmSubModelBase(film),
    ejectedMass_(0.0)
{}


ejectionModel::ejectionModel
(
    const word& modelType,
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    filmSubModelBase(film, dict, typeName, modelType),
    ejectedMass_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ejectionModel::~ejectionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void ejectionModel::correct()
{
    if (writeTime())
    {
        scalar ejectedMass0 = getModelProperty<scalar>("ejectedMass");
        ejectedMass0 += returnReduce(ejectedMass_, sumOp<scalar>());
        setModelProperty<scalar>("ejectedMass", ejectedMass0);
        ejectedMass_ = 0.0;
    }
}


scalar ejectionModel::ejectedMassTotal() const
{
    scalar ejectedMass0 = getModelProperty<scalar>("ejectedMass");
    return ejectedMass0 + returnReduce(ejectedMass_, sumOp<scalar>());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
