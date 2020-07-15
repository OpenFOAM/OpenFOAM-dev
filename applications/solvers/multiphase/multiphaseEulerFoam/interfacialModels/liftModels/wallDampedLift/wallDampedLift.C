/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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

#include "wallDampedLift.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{
    defineTypeNameAndDebug(wallDamped, 0);
    addToRunTimeSelectionTable(liftModel, wallDamped, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModels::wallDamped::wallDamped
(
    const dictionary& dict,
    const phasePair& pair
)
:
    liftModel(dict, pair),
    liftModel_(liftModel::New(dict.subDict("lift"), pair)),
    wallDampingModel_
    (
        wallDampingModel::New(dict.subDict("wallDamping"), pair)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::wallDamped::~wallDamped()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::liftModels::wallDamped::Cl() const
{
    return wallDampingModel_->damping()*liftModel_->Cl();
}


Foam::tmp<Foam::volVectorField> Foam::liftModels::wallDamped::Fi() const
{
    return wallDampingModel_->damping()*liftModel_->Fi();
}


Foam::tmp<Foam::volVectorField> Foam::liftModels::wallDamped::F() const
{
    return wallDampingModel_->damping()*liftModel_->F();
}


Foam::tmp<Foam::surfaceScalarField> Foam::liftModels::wallDamped::Ff() const
{
    return wallDampingModel_->dampingf()*liftModel_->Ff();
}


// ************************************************************************* //
