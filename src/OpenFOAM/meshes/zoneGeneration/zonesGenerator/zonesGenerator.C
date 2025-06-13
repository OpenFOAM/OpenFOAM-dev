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

#include "zonesGenerator.H"
#include "polyMesh.H"
#include "Time.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zonesGenerator, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::zonesGenerator::io(const polyMesh& mesh) const
{
    typeIOobject<IOdictionary> result
    (
        typeName,
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!result.headerOk())
    {
        result.readOpt() = IOobject::NO_READ;
    }
    else
    {
        result.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
    }

    return result;
}


void Foam::zonesGenerator::generate()
{
    zoneGeneratorList::generate();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zonesGenerator::zonesGenerator
(
    const polyMesh& mesh
)
:
    DemandDrivenMeshObject
    <
        polyMesh,
        TopoChangeableMeshObject,
        zonesGenerator
    >
    (
        io(mesh),
        mesh
    ),
    zoneGeneratorList(mesh)
{
    readHeaderOk(IOstream::ASCII, typeName);
    generate();
    addWatch();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::zonesGenerator::readData(Istream& is)
{
    is >> *this;
    zoneGeneratorList::read(*this);
    return !is.bad();
}


bool Foam::zonesGenerator::writeData(Ostream& os) const
{
    dictionary::write(os, false);
    return os.good();
}


bool Foam::zonesGenerator::read()
{
    if (regIOobject::read())
    {
        generate();
    }

    return true;
}


bool Foam::zonesGenerator::movePoints()
{
    zoneGeneratorList::movePoints();
    return true;
}


void Foam::zonesGenerator::distribute(const polyDistributionMap&)
{
    // Zones are automatically distributed
}


void Foam::zonesGenerator::topoChange(const polyTopoChangeMap&)
{
    // Regenerate zones following topology change
    generate();
}


void Foam::zonesGenerator::mapMesh(const polyMeshMap&)
{
    // Zones are automatically mapped but should be regenerated
    generate();
}


// ************************************************************************* //
