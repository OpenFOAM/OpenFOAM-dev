/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "distributor_fvMeshDistributor.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "polyDistributionMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshDistributors
{
    defineTypeNameAndDebug(distributor, 0);
    addToRunTimeSelectionTable
    (
        fvMeshDistributor,
        distributor,
        fvMesh
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fvMeshDistributors::distributor::distribute
(
    const labelList& distribution
)
{
    fvMesh& mesh = this->mesh();

    // Mesh distribution engine
    fvMeshDistribute distributor(mesh);

    // Do actual sending/receiving of mesh
    autoPtr<polyDistributionMap> map
    (
        distributor.distribute(distribution)
    );

    // Distribute the mesh data
    mesh.distribute(map);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshDistributors::distributor::distributor
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    fvMeshDistributor(mesh),
    distributor_
    (
        decompositionMethod::NewDistributor
        (
            decompositionMethod::decomposeParDict(mesh.time())
        )
    ),
    redistributionInterval_(dict.lookupOrDefault("redistributionInterval", 10)),
    maxImbalance_(dict.lookupOrDefault<scalar>("maxImbalance", 0.1)),
    timeIndex_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshDistributors::distributor::~distributor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshDistributors::distributor::update()
{
    const fvMesh& mesh = this->mesh();

    bool redistributed = false;

    if
    (
        Pstream::nProcs() > 1
     && mesh.time().timeIndex() > 1
     && timeIndex_ != mesh.time().timeIndex()
     && mesh.time().timeIndex() % redistributionInterval_ == 0
    )
    {
        timeIndex_ = mesh.time().timeIndex();

        const scalar idealNCells =
            mesh.globalData().nTotalCells()/Pstream::nProcs();

        const scalar imbalance = returnReduce
        (
            mag(1 - mesh.nCells()/idealNCells),
            maxOp<scalar>()
        );

        if (imbalance > maxImbalance_)
        {
            Info<< "Redistributing mesh with imbalance " << imbalance << endl;

            // Create new decomposition distribution
            const labelList distribution
            (
                distributor_->decompose(mesh, scalarField())
            );

            distribute(distribution);

            redistributed = true;
        }
    }

    return redistributed;
}


void Foam::fvMeshDistributors::distributor::topoChange(const polyTopoChangeMap&)
{}


void Foam::fvMeshDistributors::distributor::mapMesh(const polyMeshMap&)
{}


void Foam::fvMeshDistributors::distributor::distribute
(
    const polyDistributionMap&
)
{}


bool Foam::fvMeshDistributors::distributor::write(const bool write) const
{
    return true;
}


// ************************************************************************* //
