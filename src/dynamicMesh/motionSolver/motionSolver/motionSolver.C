/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "motionSolver.H"
#include "Time.H"
#include "polyMesh.H"
#include "dlLibraryTable.H"
#include "twoDPointCorrector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionSolver, 0);
    defineRunTimeSelectionTable(motionSolver, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::IOobject Foam::motionSolver::stealRegistration
(
    const IOdictionary& dict
)
{
    IOobject io(dict);
    if (dict.registerObject())
    {
        // De-register if necessary
        const_cast<IOdictionary&>(dict).checkOut();

        io.registerObject() = true;
    }

    return io;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionSolver::motionSolver(const polyMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "dynamicMeshDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::AUTO_WRITE
        )
    ),
    mesh_(mesh)
{}


Foam::motionSolver::motionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const word& type
)
:
    IOdictionary(stealRegistration(dict), dict),
    mesh_(mesh),
    coeffDict_(dict.subDict(type + "Coeffs"))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::motionSolver> Foam::motionSolver::New
(
    const polyMesh& mesh,
    const IOdictionary& solverDict
)
{
    const word solverTypeName(solverDict.lookup("solver"));

    Info<< "Selecting motion solver: " << solverTypeName << endl;

    const_cast<Time&>(mesh.time()).libs().open
    (
        solverDict,
        "motionSolverLibs",
        dictionaryConstructorTablePtr_
    );

    if (!dictionaryConstructorTablePtr_)
    {
        FatalErrorIn
        (
            "motionSolver::New(const polyMesh& mesh)"
        )   << "solver table is empty"
            << exit(FatalError);
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solverTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "motionSolver::New(const polyMesh&)"
        )   << "Unknown solver type "
            << solverTypeName << nl << nl
            << "Valid solver types are:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<motionSolver>(cstrIter()(mesh, solverDict));
}


Foam::autoPtr<Foam::motionSolver> Foam::motionSolver::New(const polyMesh& mesh)
{
    IOdictionary solverDict
    (
        IOobject
        (
            "dynamicMeshDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::AUTO_WRITE
        )
    );

    return New(mesh, solverDict);
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


void Foam::motionSolver::updateMesh(const mapPolyMesh& mpm)
{}


bool Foam::motionSolver::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    return true;
}


bool Foam::motionSolver::read()
{
    if (regIOobject::read())
    {
        coeffDict_ = subDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
