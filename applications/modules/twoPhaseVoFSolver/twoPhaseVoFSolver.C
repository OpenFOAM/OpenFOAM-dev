/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

#include "twoPhaseVoFSolver.H"
#include "fvcAverage.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(twoPhaseVoFSolver, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::twoPhaseVoFSolver::correctCoNum()
{
    twoPhaseSolver::correctCoNum();

    const scalarField sumPhi
    (
        interface.nearInterface()().primitiveField()
       *fvc::surfaceSum(mag(phi))().primitiveField()
    );

    alphaCoNum =
        0.5*gMax(sumPhi/mesh.V().primitiveField())*runTime.deltaTValue();

    const scalar meanAlphaCoNum =
        0.5
       *(gSum(sumPhi)/gSum(mesh.V().primitiveField()))
       *runTime.deltaTValue();

    Info<< "Interface Courant Number mean: " << meanAlphaCoNum
        << " max: " << alphaCoNum << endl;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::solvers::twoPhaseVoFSolver::correctInterface()
{
    interface.correct();
}


Foam::tmp<Foam::surfaceScalarField>
Foam::solvers::twoPhaseVoFSolver::surfaceTensionForce() const
{
    return interface.surfaceTensionForce();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::twoPhaseVoFSolver::twoPhaseVoFSolver
(
    fvMesh& mesh,
    autoPtr<twoPhaseVoFMixture> mixturePtr
)
:
    twoPhaseSolver(mesh, mixturePtr),

    interface(mixture, alpha1, alpha2, U)
{
    if (transient())
    {
        correctCoNum();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::twoPhaseVoFSolver::~twoPhaseVoFSolver()
{}


// ************************************************************************* //
