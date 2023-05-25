/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "incompressibleMultiphaseVoF.H"
#include "fvCorrectPhi.H"
#include "geometricZeroField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(incompressibleMultiphaseVoF, 0);
    addToRunTimeSelectionTable(solver, incompressibleMultiphaseVoF, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::incompressibleMultiphaseVoF::incompressibleMultiphaseVoF
(
    fvMesh& mesh
)
:
    multiphaseVoFSolver
    (
        mesh,
        autoPtr<multiphaseVoFMixture>
        (
            new incompressibleMultiphaseVoFMixture(mesh)
        )
    ),

    mixture
    (
        refCast<incompressibleMultiphaseVoFMixture>
        (
            multiphaseVoFSolver::mixture
        )
    ),

    phases(mixture.phases()),

    p
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*buoyancy.gh
    ),

    pressureReference_
    (
        p,
        p_rgh,
        pimple.dict(),
        false
    ),

    momentumTransport_
    (
        incompressible::momentumTransportModel::New
        (
            U,
            phi,
            mixture
        )
    ),

    momentumTransport(momentumTransport_())
{
    // Read the controls
    readControls();

    if (correctPhi || mesh.topoChanging())
    {
        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 1)
        );
    }

    if (!runTime.restart() || !divergent())
    {
        correctUphiBCs(U_, phi_, true);

        fv::correctPhi
        (
            phi_,
            U,
            p_rgh,
            rAU,
            autoPtr<volScalarField>(),
            pressureReference(),
            pimple
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::incompressibleMultiphaseVoF::~incompressibleMultiphaseVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleMultiphaseVoF::prePredictor()
{
    multiphaseVoFSolver::prePredictor();

    if (pimple.predictTransport())
    {
        momentumTransport.predict();
    }
}


void Foam::solvers::incompressibleMultiphaseVoF::thermophysicalPredictor()
{}


void Foam::solvers::incompressibleMultiphaseVoF::postCorrector()
{
    if (pimple.correctTransport())
    {
        momentumTransport.correct();
    }
}


// ************************************************************************* //
