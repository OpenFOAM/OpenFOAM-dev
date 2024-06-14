/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

#include "cpuLoad.H"
#include "polyDistributionMap.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "cellMapper.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cpuLoad, 0);
    optionalCpuLoad optionalCpuLoad::optionalCpuLoad_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cpuLoad::cpuLoad( const word& name, const fvMesh& mesh)
:
    DemandDrivenMeshObject<fvMesh, TopoChangeableMeshObject, cpuLoad>
    (
        name,
        mesh
    ),
    scalarField(mesh.nCells(), 0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cpuLoad::~cpuLoad()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::optionalCpuLoad& Foam::optionalCpuLoad::New
(
    const word& name,
    const fvMesh& mesh,
    const bool loadBalancing
)
{
    if (loadBalancing)
    {
        return DemandDrivenMeshObject
        <
            fvMesh,
            TopoChangeableMeshObject,
            cpuLoad
        >::New(name, mesh);
    }
    else
    {
        return optionalCpuLoad::optionalCpuLoad_;
    }
}


Foam::optionalCpuLoad& Foam::optionalCpuLoad::New
(
    const word& name,
    const polyMesh& mesh,
    const bool loadBalancing
)
{
    if (loadBalancing && isA<fvMesh>(mesh))
    {
        return New(name, refCast<const fvMesh>(mesh), loadBalancing);
    }
    else
    {
        return optionalCpuLoad::optionalCpuLoad_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cpuLoad::resetCpuTime()
{
    cpuTime_.cpuTimeIncrement();
}


void Foam::cpuLoad::cpuTimeIncrement(const label celli)
{
    operator[](celli) += cpuTime_.cpuTimeIncrement();
}


void Foam::cpuLoad::reset()
{
    scalarField::operator=(0);
}


bool Foam::cpuLoad::movePoints()
{
    return true;
}


void Foam::cpuLoad::topoChange(const polyTopoChangeMap& map)
{
    const cellMapper cellMap(map);
    cellMap(*this, *this);
}


void Foam::cpuLoad::mapMesh(const polyMeshMap& map)
{
    setSize(map.mesh().nCells());
    reset();
}


void Foam::cpuLoad::distribute(const polyDistributionMap& map)
{
    setSize(map.mesh().nCells());
    reset();
}


// ************************************************************************* //
