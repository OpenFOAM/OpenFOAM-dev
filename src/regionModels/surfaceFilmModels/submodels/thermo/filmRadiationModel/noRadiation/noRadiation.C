/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "noRadiation.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noRadiation, 0);

addToRunTimeSelectionTable
(
    filmRadiationModel,
    noRadiation,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noRadiation::noRadiation
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    filmRadiationModel(film)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noRadiation::~noRadiation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void noRadiation::correct()
{}


tmp<volScalarField::Internal> noRadiation::Shs()
{
    return volScalarField::Internal::New
    (
        IOobject::modelName("Shs", typeName),
        film().regionMesh(),
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
