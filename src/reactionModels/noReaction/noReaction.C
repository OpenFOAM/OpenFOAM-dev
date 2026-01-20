/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "noReaction.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionModels
{
    defineTypeNameAndDebug(noReaction, 0);
    addToRunTimeSelectionTable(reactionModel, noReaction, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionModels::noReaction::noReaction
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& reactionProperties
)
:
    reactionModel(modelType, thermo, turb)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionModels::noReaction::~noReaction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reactionModels::noReaction::correct()
{}


Foam::tmp<Foam::volScalarField::Internal>
Foam::reactionModels::noReaction::R(const label speciei) const
{
    return
        volScalarField::Internal::New
        (
            typedName("R_" + this->thermo().Y()[speciei].name()),
            this->mesh(),
            dimensionedScalar(dimDensity/dimTime, 0)
        );
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::reactionModels::noReaction::R(volScalarField& Y) const
{
    return tmp<fvScalarMatrix>(new fvScalarMatrix(Y, dimMass/dimTime));
}


Foam::tmp<Foam::volScalarField>
Foam::reactionModels::noReaction::Qdot() const
{
    return volScalarField::New
    (
        this->thermo().phasePropertyName(typedName("Qdot")),
        this->mesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
    );
}


bool Foam::reactionModels::noReaction::read()
{
    if (reactionModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
