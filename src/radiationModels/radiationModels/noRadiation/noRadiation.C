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
#include "physicoChemicalConstants.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
    defineTypeNameAndDebug(noRadiation, 0);
    addToRadiationRunTimeSelectionTables(noRadiation);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::noRadiation::noRadiation(const volScalarField& T)
:
    radiationModel(T)
{}


Foam::radiationModels::noRadiation::noRadiation
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(T)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::noRadiation::~noRadiation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiationModels::noRadiation::correct()
{}


void Foam::radiationModels::noRadiation::calculate()
{}


bool Foam::radiationModels::noRadiation::read()
{
    return radiationModel::read();
}


Foam::tmp<Foam::volScalarField> Foam::radiationModels::noRadiation::Rp() const
{
    return volScalarField::New
    (
        "Rp",
        mesh_,
        dimensionedScalar
        (
            constant::physicoChemical::sigma.dimensions()/dimLength,
            0
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiationModels::noRadiation::Ru() const
{
    return volScalarField::Internal::New
    (
        "Ru",
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    );
}


// ************************************************************************* //
