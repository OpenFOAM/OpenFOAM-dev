/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "cloud_functionObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloud, 0);
}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloud::~cloud()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::cloud::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::cloud::executeAtStart() const
{
    return false;
}


bool Foam::functionObjects::cloud::execute()
{
    cloudPtr_->solve();
    Info<< endl;
    return true;
}


bool Foam::functionObjects::cloud::write()
{
    cloudPtr_->mesh().write();
    return true;
}


void Foam::functionObjects::cloud::movePoints(const polyMesh&)
{}


void Foam::functionObjects::cloud::topoChange(const polyTopoChangeMap& map)
{
    cloudPtr_->topoChange(map);
}


void Foam::functionObjects::cloud::mapMesh(const polyMeshMap& map)
{
    cloudPtr_->mapMesh(map);
}


void Foam::functionObjects::cloud::distribute(const polyDistributionMap& map)
{
    cloudPtr_->distribute(map);
}


// ************************************************************************* //
