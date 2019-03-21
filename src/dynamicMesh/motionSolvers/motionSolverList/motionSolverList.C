/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "motionSolverList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionSolverList, 0);
    addToRunTimeSelectionTable
    (
        motionSolver,
        motionSolverList,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionSolverList::motionSolverList
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    motionSolver(mesh, dict, typeName),
    motionSolvers_
    (
        dict.lookup("solvers"),
        motionSolver::iNew(mesh)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionSolverList::~motionSolverList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::motionSolverList::curPoints() const
{
    if (motionSolvers_.size())
    {
        // Accumulated displacement
        pointField disp(motionSolvers_[0].curPoints() - mesh().points());

        for (label i = 1; i < motionSolvers_.size(); i++)
        {
            disp += motionSolvers_[i].curPoints() - mesh().points();
        }

        return mesh().points() + disp;
    }
    else
    {
        return mesh().points();
    }
}


void Foam::motionSolverList::solve()
{
    forAll(motionSolvers_, i)
    {
        motionSolvers_[i].solve();
    }
}


void Foam::motionSolverList::movePoints(const pointField& points)
{
    forAll(motionSolvers_, i)
    {
        motionSolvers_[i].movePoints(points);
    }
}


void Foam::motionSolverList::updateMesh(const mapPolyMesh& mpm)
{
    forAll(motionSolvers_, i)
    {
        motionSolvers_[i].updateMesh(mpm);
    }
}


// ************************************************************************* //
