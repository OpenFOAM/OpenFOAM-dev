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

#include "fvMeshDistributorsLoadBalancer.H"
#include "decompositionMethod.H"
#include "cpuLoad.H"
#include "globalMeshData.H"
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
     && mesh.time().timeIndex() - mesh.time().startTimeIndex() > 1
     && timeIndex_ != mesh.time().timeIndex()
    )
    {
        timeIndex_ = mesh.time().timeIndex();

        // Get the CPU time fer this processor which includes waiting time
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

            scalarList procCpuLoads(cpuLoads.size());

            label l = 0;
            forAllConstIter(HashTable<cpuLoad*>, cpuLoads, iter)
            {
                procCpuLoads[l++] = sum(*iter());
            }

            List<scalarList> allProcCpuLoads(Pstream::nProcs());
            allProcCpuLoads[Pstream::myProcNo()] = procCpuLoads;
            Pstream::gatherList(allProcCpuLoads);
            Pstream::scatterList(allProcCpuLoads);

            scalarList sumProcCpuLoads(procCpuLoads.size(), scalar(0));
            scalarList maxProcCpuLoads(procCpuLoads.size(), scalar(0));
            forAll(maxProcCpuLoads, l)
            {
                forAll(allProcCpuLoads, proci)
                {
                    sumProcCpuLoads[l] += allProcCpuLoads[proci][l];

                    maxProcCpuLoads[l] =
                        max(maxProcCpuLoads[l], allProcCpuLoads[proci][l]);
                }
            }

            // Sum over loads of the maximum load CPU time per processor
            const scalar sumMaxProcCpuLoad(sum(maxProcCpuLoads));

            // Maximum number of cells per processor
            const label maxNcells = returnReduce(mesh.nCells(), maxOp<label>());

            // Maximum processor CPU time spent doing basic CFD
            const scalar maxBaseCpuTime =
                returnReduce(timeStepCpuTime, maxOp<scalar>())
              - sumMaxProcCpuLoad;

            const scalar cellBaseCpuTime = maxBaseCpuTime/maxNcells;

            // Processor CPU time spent doing basic CFD, not waiting
            const scalar baseCpuTime = mesh.nCells()*cellBaseCpuTime;

            // Maximum total CPU time
            const scalar maxProcCpuTime = maxBaseCpuTime + sumMaxProcCpuLoad;

            // Total CPU time for this processor not waiting
            const scalar procCpuTime = baseCpuTime + sum(procCpuLoads);

            // Average processor CPU time
            const scalar averageProcessorCpuTime =
                returnReduce(procCpuTime, sumOp<scalar>())/Pstream::nProcs();

            const scalar imbalance =
                (maxProcCpuTime - averageProcessorCpuTime)
               /averageProcessorCpuTime;

            Info<< nl << type() << nl;

            l = 0;
            forAllConstIter(HashTable<cpuLoad*>, cpuLoads, iter)
            {
                Info<< "    Imbalance of load " << iter()->name() << ": "
                    << (
                          maxProcCpuLoads[l]
                        - sumProcCpuLoads[l]/Pstream::nProcs()
                       )/averageProcessorCpuTime
                    << endl;
                l++;
            }

            Info<< "    Imbalance of base load " << ": "
                << (
                      maxBaseCpuTime
                    - mesh.globalData().nTotalCells()*cellBaseCpuTime
                     /Pstream::nProcs()
                   )/averageProcessorCpuTime
                << endl;

            Info<< "    Total imbalance " << imbalance << endl;

            if (imbalance > maxImbalance_)
            {
                Info<< "    Redistributing mesh" << endl;

                scalarField weights;

                if (multiConstraint_)
                {
                    const label nWeights = cpuLoads.size() + 1;

                    weights.setSize(nWeights*mesh.nCells());

                    for (label i=0; i<mesh.nCells(); i++)
                    {
                        weights[nWeights*i] = cellBaseCpuTime;
                    }

                    label l = 1;
                    forAllConstIter(HashTable<cpuLoad*>, cpuLoads, iter)
                    {
                        const scalarField& cpuLoadField = *iter();

                        forAll(cpuLoadField, i)
                        {
                            weights[nWeights*i + l] = cpuLoadField[i];
                        }

                        iter()->checkOut();

                        l++;
                    }
                }
                else
                {
                    weights.setSize(mesh.nCells(), cellBaseCpuTime);

                    forAllConstIter(HashTable<cpuLoad*>, cpuLoads, iter)
                    {
                        weights += *iter();
                        iter()->checkOut();
                    }
                }

                // Create new decomposition distribution
                const labelList distribution
                (
                    distributor_->decompose(mesh, weights)
                );

                distribute(distribution);

                redistributed = true;

                Info<< endl;
            }
            else
            {
                forAllIter(HashTable<cpuLoad*>, cpuLoads, iter)
                {
                    iter()->checkOut();
                }
            }
        }
        else
        {
            forAllIter(HashTable<cpuLoad*>, cpuLoads, iter)
            {
                iter()->checkOut();
            }
        }
    }

    return redistributed;
}


// ************************************************************************* //
