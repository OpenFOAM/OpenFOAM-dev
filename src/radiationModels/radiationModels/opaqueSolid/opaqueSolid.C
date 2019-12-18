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

#include "opaqueSolid.H"
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
    defineTypeNameAndDebug(opaqueSolid, 0);

    addToRadiationRunTimeSelectionTables(opaqueSolid);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::opaqueSolid::opaqueSolid(const volScalarField& T)
:
    radiationModel(typeName, T)
{}


Foam::radiationModels::opaqueSolid::opaqueSolid
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::opaqueSolid::~opaqueSolid()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiationModels::opaqueSolid::read()
{
    return radiationModel::read();
}


void Foam::radiationModels::opaqueSolid::calculate()
{}


Foam::tmp<Foam::volScalarField> Foam::radiationModels::opaqueSolid::Rp() const
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
Foam::radiationModels::opaqueSolid::Ru() const
{
    return volScalarField::Internal::New
    (
        "Ru",
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    );
}


// ************************************************************************* //
