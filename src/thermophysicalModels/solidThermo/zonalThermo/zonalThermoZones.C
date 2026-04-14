/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "zonalThermoZones.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zonalThermoZones, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::zonalThermoZones::update()
{
    cellZones_.resize(mesh().nCells());
    cellZones_ = -1;

    forAll(zones_, zonei)
    {
        const cellZone& cz = mesh().cellZones()[zones_[zonei]];

        forAll(cz, i)
        {
            if (cellZones_[cz[i]] != -1)
            {
                FatalErrorInFunction
                    << "Cell zones must not overlap"
                    << exit(FatalError);
            }

            cellZones_[cz[i]] = zonei;
        }
    }

    forAll(cellZones_, celli)
    {
        if (cellZones_[celli] == -1)
        {
            FatalErrorInFunction
                << "Cell zones must span the entire mesh"
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zonalThermoZones::zonalThermoZones
(
    const word& name,
    const fvMesh& mesh,
    const hashedWordList& zones
)
:
    DemandDrivenMeshObject<fvMesh, TopoChangeableMeshObject, zonalThermoZones>
    (
        name,
        mesh
    ),
    zones_(zones),
    cellZones_()
{
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zonalThermoZones::~zonalThermoZones()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::zonalThermoZones::movePoints()
{
    return true;
}


void Foam::zonalThermoZones::topoChange(const polyTopoChangeMap&)
{
    update();
}


void Foam::zonalThermoZones::mapMesh(const polyMeshMap&)
{
    update();
}


void Foam::zonalThermoZones::distribute(const polyDistributionMap&)
{
    update();
}


// ************************************************************************* //
