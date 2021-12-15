/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "fvMeshDistributorsDistributor.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"
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

void Foam::fvMeshDistributors::distributor::readDict()
{
    const dictionary& distributorDict(dict());

    redistributionInterval_ =
        distributorDict.lookupOrDefault("redistributionInterval", 10);

    maxImbalance_ =
        distributorDict.lookupOrDefault<scalar>("maxImbalance", 0.1);
}


void Foam::fvMeshDistributors::distributor::distribute()
{
    fvMesh& mesh = this->mesh();

    // Create new decomposition distribution
    labelList distribution
    (
        distributor_->decompose(mesh, mesh.cellCentres())
    );

    // Mesh distribution engine
    fvMeshDistribute distributor(mesh);

    // Do actual sending/receiving of mesh
    autoPtr<mapDistributePolyMesh> map
    (
        distributor.distribute(distribution)
    );

    // Distribute the mesh data
    mesh.distribute(map);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshDistributors::distributor::distributor(fvMesh& mesh)
:
    fvMeshDistributor(mesh),
    distributor_
    (
        decompositionMethod::NewDistributor
        (
            decompositionMethod::decomposeParDict(mesh.time())
        )
    ),
    redistributionInterval_(1),
    maxImbalance_(0.1),
    timeIndex_(-1)
{
    readDict();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshDistributors::distributor::~distributor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshDistributors::distributor::update()
{
    if
    (
        Pstream::nProcs() > 1
     && mesh().time().timeIndex() > 1
     && timeIndex_ != mesh().time().timeIndex()
     && mesh().time().timeIndex() % redistributionInterval_ == 0
    )
    {
        timeIndex_ = mesh().time().timeIndex();

        const scalar idealNCells =
            mesh().globalData().nTotalCells()/Pstream::nProcs();

        const scalar imbalance = returnReduce
        (
            mag(1 - mesh().nCells()/idealNCells),
            maxOp<scalar>()
        );

        if (imbalance > maxImbalance_)
        {
            Info<< "Redistributing mesh with imbalance " << imbalance << endl;

            distribute();

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


void Foam::fvMeshDistributors::distributor::updateMesh(const mapPolyMesh&)
{}


void Foam::fvMeshDistributors::distributor::distribute
(
    const mapDistributePolyMesh&
)
{}


bool Foam::fvMeshDistributors::distributor::write(const bool write) const
{
    return true;
}


// ************************************************************************* //
