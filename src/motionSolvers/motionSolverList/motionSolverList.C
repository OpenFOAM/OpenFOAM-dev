/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2024 OpenFOAM Foundation
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
#include "polyMesh.H"
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
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    motionSolver(name, mesh, typeName)
{
    const dictionary& solversDict = dict.subDict("solvers");

    forAllConstIter(dictionary, solversDict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& dict = iter().dict();

            motionSolvers_.insert
            (
                name,
                motionSolver::New(name, mesh, dict).ptr()
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionSolverList::~motionSolverList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::motionSolverList::curPoints() const
{
    if (motionSolvers_.size())
    {
        // Accumulated displacement
        pointField disp(mesh().nPoints(), Zero);

        forAllConstIter(PtrDictionary<motionSolver>, motionSolvers_, iter)
        {
            disp += iter().curPoints() - mesh().points();
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
    forAllIter(PtrDictionary<motionSolver>, motionSolvers_, iter)
    {
        iter().solve();
    }
}


void Foam::motionSolverList::topoChange(const polyTopoChangeMap& map)
{
    forAllIter(PtrDictionary<motionSolver>, motionSolvers_, iter)
    {
        iter().topoChange(map);
    }
}


void Foam::motionSolverList::mapMesh(const polyMeshMap& map)
{
    forAllIter(PtrDictionary<motionSolver>, motionSolvers_, iter)
    {
        iter().mapMesh(map);
    }
}


void Foam::motionSolverList::distribute
(
    const polyDistributionMap& map
)
{
    forAllIter(PtrDictionary<motionSolver>, motionSolvers_, iter)
    {
        iter().distribute(map);
    }
}


void Foam::motionSolverList::movePoints(const pointField& points)
{
    forAllIter(PtrDictionary<motionSolver>, motionSolvers_, iter)
    {
        iter().movePoints(points);
    }
}


// ************************************************************************* //
