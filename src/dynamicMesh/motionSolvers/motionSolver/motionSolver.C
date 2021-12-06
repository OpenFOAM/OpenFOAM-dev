/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
#include "twoDPointCorrector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionSolver, 0);
    defineRunTimeSelectionTable(motionSolver, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionSolver::motionSolver
(
    const polyMesh& mesh,
    const dictionary& dict,
    const word& type
)
:
    mesh_(mesh),
    coeffDict_(dict.optionalSubDict(type + "Coeffs"))
{}


Foam::autoPtr<Foam::motionSolver> Foam::motionSolver::clone() const
{
    NotImplemented;
    return autoPtr<motionSolver>(nullptr);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::motionSolver> Foam::motionSolver::New
(
    const polyMesh& mesh,
    const dictionary& solverDict
)
{
    if (solverDict.found("motionSolvers"))
    {
        return autoPtr<motionSolver>(new motionSolverList(mesh, solverDict));
    }
    else
    {
        const word solverTypeName = solverDict.lookup<word>("motionSolver");

        Info<< "Selecting motion solver: " << solverTypeName << endl;

        libs.open
        (
            solverDict,
            "motionSolverLibs",
            dictionaryConstructorTablePtr_
        );

        if (!dictionaryConstructorTablePtr_)
        {
            FatalErrorInFunction
                << "solver table is empty"
                << exit(FatalError);
        }

        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(solverTypeName);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown solver type "
                << solverTypeName << nl << nl
                << "Valid solver types are:" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<motionSolver>(cstrIter()(mesh, solverDict));
    }
}


Foam::motionSolver::iNew::iNew(const polyMesh& mesh)
:
    mesh_(mesh)
{}


Foam::autoPtr<Foam::motionSolver> Foam::motionSolver::iNew::operator()
(
    Istream& is
) const
{
    dictionaryEntry dict(dictionary::null, is);

    return motionSolver::New(mesh_, dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionSolver::~motionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::motionSolver::newPoints()
{
    solve();
    return curPoints();
}


void Foam::motionSolver::twoDCorrectPoints(pointField& p) const
{
    twoDPointCorrector::New(mesh_).correctPoints(p);
}


bool Foam::motionSolver::write() const
{
    return true;
}


// ************************************************************************* //
