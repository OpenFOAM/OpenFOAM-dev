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

Foam::cpuLoad::cpuLoad(const fvMesh& mesh, const word& name)
:
    DemandDrivenMeshObject<fvMesh, TopoChangeableMeshObject, cpuLoad>
    (
        mesh,
        IOobject
        (
            name,
            mesh.time().name(),
            mesh
        )
    ),
    scalarField(mesh.nCells(), 0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cpuLoad::~cpuLoad()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::optionalCpuLoad& Foam::optionalCpuLoad::New
(
    const fvMesh& mesh,
    const word& name,
    const bool loadBalancing
)
{
    if (loadBalancing)
    {
        if
        (
            mesh.thisDb().objectRegistry::template
            foundObject<cpuLoad>
            (
                name
            )
        )
        {
            return mesh.thisDb().objectRegistry::template
            lookupObjectRef<cpuLoad>
            (
                name
            );
        }
        else
        {
            if (cpuLoad::debug)
            {
                InfoInFunction
                    << "constructing " << name
                    << " for region " << mesh.name() << endl;
            }

            cpuLoad* cpuLoadPtr(new cpuLoad(mesh, name));

            regIOobject::store(cpuLoadPtr);

            return *cpuLoadPtr;
        }
    }
    else
    {
        return optionalCpuLoad::optionalCpuLoad_;
    }
}


Foam::optionalCpuLoad& Foam::optionalCpuLoad::New
(
    const polyMesh& mesh,
    const word& name,
    const bool loadBalancing
)
{
    if (loadBalancing && isA<fvMesh>(mesh))
    {
        return New(refCast<const fvMesh>(mesh), name, loadBalancing);
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
