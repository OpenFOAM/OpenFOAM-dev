/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "noPhaseChange.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmSubModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noPhaseChange, 0);
addToRunTimeSelectionTable(phaseChangeModel, noPhaseChange, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noPhaseChange::noPhaseChange
(
    surfaceFilm& film,
    const dictionary&
)
:
    phaseChangeModel(film)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noPhaseChange::~noPhaseChange()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void noPhaseChange::correct
(
    const scalar dt,
    scalarField& availableMass,
    volScalarField& dMass,
    volScalarField& dEnergy
)
{}


void noPhaseChange::correctModel
(
    const scalar,
    scalarField&,
    scalarField&,
    scalarField&
)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmSubModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
