/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "multiphaseEuler.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::solvers::multiphaseEuler::cellMomentumPredictor()
{
    Info<< "Constructing momentum equations" << endl;

    phaseSystem& fluid(fluid_);

    autoPtr<phaseSystem::momentumTransferTable>
        momentumTransferPtr(fluid.momentumTransfer());

    phaseSystem::momentumTransferTable&
        momentumTransfer(momentumTransferPtr());

    const PtrList<volScalarField> Kds(fluid.Kds());

    forAll(fluid.movingPhases(), movingPhasei)
    {
        phaseModel& phase = fluid.movingPhases()[movingPhasei];

        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();
        volVectorField& U = phase.URef();

        UEqns.set
        (
            phase.index(),
            new fvVectorMatrix
            (
                phase.UEqn()
             ==
               *momentumTransfer[phase.name()]
              + fvModels().source(alpha, rho, U)
                // - fvm::Sp(Kds[phase.index()], U)
            )
        );

        UEqns[phase.index()].relax();
        fvConstraints().constrain(UEqns[phase.index()]);
        U.correctBoundaryConditions();
        fvConstraints().constrain(U);
    }
}


void Foam::solvers::multiphaseEuler::faceMomentumPredictor()
{
    Info<< "Constructing face momentum equations" << endl;

    phaseSystem& fluid(fluid_);

    autoPtr<phaseSystem::momentumTransferTable>
        momentumTransferPtr(fluid.momentumTransferf());

    phaseSystem::momentumTransferTable&
        momentumTransfer(momentumTransferPtr());

    forAll(fluid.movingPhases(), movingPhasei)
    {
        phaseModel& phase = fluid.movingPhases()[movingPhasei];

        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();
        volVectorField& U = phase.URef();

        UEqns.set
        (
            phase.index(),
            new fvVectorMatrix
            (
                phase.UfEqn()
             ==
               *momentumTransfer[phase.name()]
              + fvModels().source(alpha, rho, U)
            )
        );

        UEqns[phase.index()].relax();
        fvConstraints().constrain(UEqns[phase.index()]);
        U.correctBoundaryConditions();
        fvConstraints().constrain(U);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::multiphaseEuler::momentumPredictor()
{
    if (pimple.flow())
    {
        UEqns.setSize(phases.size());

        if (faceMomentum)
        {
            faceMomentumPredictor();
        }
        else
        {
            cellMomentumPredictor();
        }
    }
}


// ************************************************************************* //
