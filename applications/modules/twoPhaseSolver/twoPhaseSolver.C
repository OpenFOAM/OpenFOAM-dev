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

#include "twoPhaseSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(twoPhaseSolver, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::twoPhaseSolver::twoPhaseSolver
(
    fvMesh& mesh,
    autoPtr<twoPhaseVoFMixture> mixturePtr
)
:
    VoFSolver(mesh, autoPtr<VoFMixture>(mixturePtr.ptr())),

    mixture(refCast<twoPhaseVoFMixture>(VoFSolver::mixture_)),

    alpha1(mixture.alpha1()),
    alpha2(mixture.alpha2()),

    alphaRestart
    (
        typeIOobject<surfaceScalarField>
        (
            IOobject::groupName("alphaPhi", alpha1.group()),
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ).headerOk()
    ),

    alphaPhi1
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", alpha1.group()),
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        phi*fvc::interpolate(alpha1)
    ),
    alphaPhi2
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", alpha2.group()),
            runTime.name(),
            mesh
        ),
        phi - alphaPhi1
    )
{
    mesh.schemes().setFluxRequired(alpha1.name());

    if (alphaRestart)
    {
        Info << "Restarting alpha" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::twoPhaseSolver::~twoPhaseSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::twoPhaseSolver::preSolve()
{
    VoFSolver::preSolve();

    // Do not apply previous time-step mesh compression flux
    // if the mesh topology changed
    if (mesh().topoChanged())
    {
        talphaPhi1Corr0.clear();
    }
}


void Foam::solvers::twoPhaseSolver::prePredictor()
{
    VoFSolver::prePredictor();
    alphaPredictor();
    mixture.correct();
}


// ************************************************************************* //
