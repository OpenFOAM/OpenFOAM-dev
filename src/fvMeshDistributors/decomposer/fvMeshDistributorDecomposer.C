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

#include "fvMeshDistributorDecomposer.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshDistributors
{
    defineTypeNameAndDebug(decomposer, 0);
    addToRunTimeSelectionTable
    (
        fvMeshDistributor,
        decomposer,
        fvMesh
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fvMeshDistributors::decomposer::readDict()
{
    const dictionary& decomposerDict(dict());

    redistributionInterval_ =
        decomposerDict.lookupOrDefault("redistributionInterval", 10);

    maxImbalance_ =
        decomposerDict.lookupOrDefault<scalar>("maxImbalance", 0.1);
}


void Foam::fvMeshDistributors::decomposer::redecompose()
{
    fvMesh& mesh = this->mesh();

    IOdictionary decompositionDict
    (
        IOobject
        (
            "decomposeParDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    labelList finalDecomp;

    // Create decompositionMethod and new decomposition
    {
        autoPtr<decompositionMethod> decomposer
        (
            decompositionMethod::New
            (
                decompositionDict
            )
        );

        if (!decomposer().parallelAware())
        {
            WarningInFunction
                << "You have selected decomposition method "
                << decomposer().typeName
                << " which does" << endl
                << "not synchronise the decomposition across"
                << " processor patches." << endl
                << "    You might want to select a decomposition method which"
                << " is aware of this. Continuing."
                << endl;
        }

        finalDecomp = decomposer().decompose(mesh, mesh.cellCentres());
    }

    // Mesh distribution engine
    fvMeshDistribute distributor(mesh);

    // Do actual sending/receiving of mesh
    autoPtr<mapDistributePolyMesh> map = distributor.distribute(finalDecomp);

    mesh.distribute(map);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshDistributors::decomposer::decomposer(fvMesh& mesh)
:
    fvMeshDistributor(mesh),
    redistributionInterval_(1),
    maxImbalance_(0.1),
    timeIndex_(-1)
{
    readDict();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshDistributors::decomposer::~decomposer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshDistributors::decomposer::update()
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
            if (debug)
            {
                Info<< "Redistributing mesh with imbalance "
                    << imbalance << endl;
            }

            redecompose();

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


void Foam::fvMeshDistributors::decomposer::updateMesh(const mapPolyMesh&)
{}


void Foam::fvMeshDistributors::decomposer::distribute
(
    const mapDistributePolyMesh&
)
{}


bool Foam::fvMeshDistributors::decomposer::write(const bool write) const
{
    return true;
}


// ************************************************************************* //
