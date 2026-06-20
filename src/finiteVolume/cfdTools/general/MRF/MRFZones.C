/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2026 OpenFOAM Foundation
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

#include "MRFZones.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MRFZones, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::MRFZones::createIOobject
(
    const fvMesh& mesh
) const
{
    typeIOobject<IOdictionary> io
    (
        "MRFProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        Info<< indentOrNl
            << "Constructing MRF zones from " << io.name()
            << endl;

        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFZones::MRFZones
(
    const fvMesh& mesh
)
:
    DemandDrivenMeshObject
    <
        fvMesh,
        TopoChangeableMeshObject,
        MRFZones,
        IOdictionary
    >
    (
        createIOobject(mesh),
        mesh
    ),
    MRFZoneList(mesh, *this)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MRFZones::preUpdateMesh()
{}


bool Foam::MRFZones::movePoints()
{
    update();
    return true;
}


void Foam::MRFZones::topoChange(const polyTopoChangeMap& map)
{
    update();
}


void Foam::MRFZones::mapMesh(const polyMeshMap& map)
{
    update();
}


void Foam::MRFZones::distribute(const polyDistributionMap& map)
{
    update();
}


bool Foam::MRFZones::read()
{
    if (regIOobject::read())
    {
        MRFZoneList::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
