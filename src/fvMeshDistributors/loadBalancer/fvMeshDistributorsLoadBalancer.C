/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "fvMeshDistributorsLoadBalancer.H"
#include "decompositionMethod.H"
#include "cpuLoad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshDistributors
{
    defineTypeNameAndDebug(loadBalancer, 0);
    addToRunTimeSelectionTable
    (
        fvMeshDistributor,
        loadBalancer,
        fvMesh
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fvMeshDistributors::loadBalancer::readDict()
{
    distributor::readDict();

    const dictionary& distributorDict(dict());

    multiConstraint_ =
        distributorDict.lookupOrDefault<Switch>("multiConstraint", true);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshDistributors::loadBalancer::loadBalancer(fvMesh& mesh)
:
    distributor(mesh)
{
    readDict();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshDistributors::loadBalancer::~loadBalancer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshDistributors::loadBalancer::update()
{
    const fvMesh& mesh = this->mesh();

    bool redistributed = false;

    if
    (
        Pstream::nProcs() > 1
     && mesh.time().timeIndex() > 1
     && timeIndex_ != mesh.time().timeIndex()
    )
    {
        timeIndex_ = mesh.time().timeIndex();

        const scalar timeStepCpuTime = cpuTime_.cpuTimeIncrement();

        // CPU loads per cell
        HashTable<cpuLoad*> cpuLoads(this->mesh().lookupClass<cpuLoad>());

        if (!cpuLoads.size())
        {
            FatalErrorInFunction
                << "No CPU loads have been allocated"
                << exit(FatalError);
        }

        if (mesh.time().timeIndex() % redistributionInterval_ == 0)
        {
            timeIndex_ = mesh.time().timeIndex();

            scalar sumCpuLoad = 0;

            forAllConstIter(HashTable<cpuLoad*>, cpuLoads, iter)
            {
                sumCpuLoad += sum(iter()->field());
            }

            const scalar cellCFDCpuTime = returnReduce
            (
                (timeStepCpuTime - sumCpuLoad)/mesh.nCells(),
                minOp<scalar>()
            );

            // Total CPU time for this processor
            const scalar processorCpuTime =
                mesh.nCells()*cellCFDCpuTime + sumCpuLoad;

            // Average processor CPU time
            const scalar averageProcessorCpuTime =
                returnReduce(processorCpuTime, sumOp<scalar>())
               /Pstream::nProcs();

            Pout<< "imbalance "
                << " " << sumCpuLoad
                << " " << mesh.nCells()*cellCFDCpuTime
                << " " << processorCpuTime
                << " " << averageProcessorCpuTime << endl;

            const scalar imbalance = returnReduce
            (
                mag(1 - processorCpuTime/averageProcessorCpuTime),
                maxOp<scalar>()
            );

            scalarField weights;

            if (multiConstraint_)
            {
                const int nWeights = cpuLoads.size() + 1;

                weights.setSize(nWeights*mesh.nCells());

                for (label i=0; i<mesh.nCells(); i++)
                {
                    weights[nWeights*i] = cellCFDCpuTime;
                }

                label loadi = 1;
                forAllConstIter(HashTable<cpuLoad*>, cpuLoads, iter)
                {
                    const scalarField& cpuLoadField = iter()->field();

                    forAll(cpuLoadField, i)
                    {
                        weights[nWeights*i + loadi] = cpuLoadField[i];
                    }

                    loadi++;
                }
            }
            else
            {
                weights.setSize(mesh.nCells(), cellCFDCpuTime);

                forAllConstIter(HashTable<cpuLoad*>, cpuLoads, iter)
                {
                    weights += iter()->field();
                }
            }

            if (imbalance > maxImbalance_)
            {
                Info<< "Redistributing mesh with imbalance "
                    << imbalance << endl;

                // Create new decomposition distribution
                const labelList distribution
                (
                    distributor_->decompose(mesh, weights)
                );

                distribute(distribution);

                redistributed = true;
            }
        }

        forAllIter(HashTable<cpuLoad*>, cpuLoads, iter)
        {
            iter()->checkOut();
        }
    }

    return redistributed;
}


// ************************************************************************* //
