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
#include "volFields.H"
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

        // Chemistry CPU load per cell
        volScalarField::Internal& chemistryCpuTimeReg =
            mesh.lookupObjectRef<volScalarField::Internal>
            (
                "chemistryCpuTime"
            );

        const scalarField& chemistryCpuTime = chemistryCpuTimeReg.field();

        if (mesh.time().timeIndex() % redistributionInterval_ == 0)
        {
            timeIndex_ = mesh.time().timeIndex();

            const scalar sumCpuLoad(sum(chemistryCpuTime));

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
                const int nWeights = 2;

                weights.setSize(nWeights*mesh.nCells());

                forAll(chemistryCpuTime, i)
                {
                    weights[nWeights*i] = cellCFDCpuTime;
                    weights[nWeights*i + 1] = chemistryCpuTime[i];
                }
            }
            else
            {
                weights = chemistryCpuTime + cellCFDCpuTime;
            }

            if (imbalance > maxImbalance_)
            {
                Info<< "Redistributing mesh with imbalance "
                    << imbalance << endl;

                // Create new decomposition distribution
                const labelList distribution
                (
                    distributor_->decompose(mesh, mesh.cellCentres(), weights)
                );

                distribute(distribution);

                redistributed = true;
            }
        }

        chemistryCpuTimeReg.checkOut();
    }

    return redistributed;
}


// ************************************************************************* //
