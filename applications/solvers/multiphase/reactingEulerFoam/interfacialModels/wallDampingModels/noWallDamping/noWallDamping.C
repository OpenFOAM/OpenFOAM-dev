/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "noWallDamping.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallDampingModels
{
    defineTypeNameAndDebug(noWallDamping, 0);
    addToRunTimeSelectionTable
    (
        wallDampingModel,
        noWallDamping,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallDampingModels::noWallDamping::noWallDamping
(
    const dictionary& dict,
    const phasePair& pair
)
:
    wallDampingModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallDampingModels::noWallDamping::~noWallDamping()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::wallDampingModels::noWallDamping::damp
(
    const tmp<volScalarField>& Cl
) const
{
    return Cl;
}


Foam::tmp<Foam::volVectorField>
Foam::wallDampingModels::noWallDamping::damp
(
    const tmp<volVectorField>& F
) const
{
    return F;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::wallDampingModels::noWallDamping::damp
(
    const tmp<surfaceScalarField>& Ff
) const
{
    return Ff;
}


// ************************************************************************* //
