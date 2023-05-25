/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ejectionModel, 0);
    defineRunTimeSelectionTable(ejectionModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ejectionModel::ejectionModel
(
    const dictionary& dict,
    const solvers::isothermalFilm& film
)
:
    film_(film),
    rate_
    (
        volScalarField::Internal::New
        (
            "ejectionRate",
            film.mesh,
            dimensionedScalar(dimless/dimTime, 0)
        )
    ),
    diameter_
    (
        volScalarField::Internal::New
        (
            "ejectionDiameter",
            film.mesh,
            dimensionedScalar(dimLength, 0)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ejectionModel::~ejectionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ejectionModel::topoChange(const polyTopoChangeMap&)
{
    // Resize in case of mesh change
    rate_.setSize(film_.mesh.nCells());
    diameter_.setSize(film_.mesh.nCells());
}


void Foam::ejectionModel::mapMesh(const polyMeshMap&)
{
    // Resize in case of mesh change
    rate_.setSize(film_.mesh.nCells());
    diameter_.setSize(film_.mesh.nCells());
}


void Foam::ejectionModel::distribute(const polyDistributionMap&)
{
    // Resize in case of mesh change
    rate_.setSize(film_.mesh.nCells());
    diameter_.setSize(film_.mesh.nCells());
}


bool Foam::ejectionModel::movePoints()
{
    return true;
}


// ************************************************************************* //
