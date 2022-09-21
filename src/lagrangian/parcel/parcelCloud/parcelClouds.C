/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "parcelClouds.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parcelClouds::parcelClouds
(
    const fvMesh& mesh,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g
)
:
    MeshObject<fvMesh, UpdateableMeshObject, parcelClouds>(mesh),
    parcelCloudList(rho, U, mu, g)
{}


Foam::parcelClouds::parcelClouds
(
    const fvMesh& mesh,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const fluidThermo& carrierThermo
)
:
    MeshObject<fvMesh, UpdateableMeshObject, parcelClouds>(mesh),
    parcelCloudList(rho, U, g, carrierThermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::parcelClouds::~parcelClouds()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::parcelClouds::preUpdateMesh()
{
    parcelCloudList::storeGlobalPositions();
}


bool Foam::parcelClouds::movePoints()
{
    return true;
}


void Foam::parcelClouds::topoChange(const polyTopoChangeMap& map)
{
    parcelCloudList::topoChange(map);
}


void Foam::parcelClouds::mapMesh(const polyMeshMap& map)
{
    parcelCloudList::mapMesh(map);
}


void Foam::parcelClouds::distribute(const polyDistributionMap& map)
{
    parcelCloudList::distribute(map);
}


// ************************************************************************* //
