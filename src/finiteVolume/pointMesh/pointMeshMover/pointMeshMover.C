/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "pointMeshMover.H"
#include "polyMesh.H"
#include "dictionaryEntry.H"
#include "twoDPointCorrector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointMeshMover, 0);
    defineRunTimeSelectionTable(pointMeshMover, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMeshMover::pointMeshMover(const polyMesh& mesh, const word& type)
:
    mesh_(mesh)
{}


Foam::autoPtr<Foam::pointMeshMover> Foam::pointMeshMover::clone() const
{
    NotImplemented;
    return autoPtr<pointMeshMover>(nullptr);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::pointMeshMover> Foam::pointMeshMover::New
(
    const polyMesh& mesh,
    const dictionary& dict
)
{
    word type = dict.lookupOrDefaultBackwardsCompatible<word>
    (
        {pointMeshMover::typeName, "motionSolver"}, word::null
    );

    if (type == word::null)
    {
        type = dict.lookup<word>("type");
    }

    Info<< indentOrNl << "Selecting pointMeshMover: " << type << endl;

    libs.open
    (
        dict,
        "libs",
        dictionaryConstructorTablePtr_
    );

    if (!dictionaryConstructorTablePtr_)
    {
        FatalIOErrorInFunction(dict)
            << "solver table is empty"
            << exit(FatalIOError);
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown solver type "
            << type << nl << nl
            << "Valid solver types are:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<pointMeshMover>
    (
        cstrIter()(mesh, dict.optionalTypeDict(type))
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMeshMover::~pointMeshMover()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointMeshMover::twoDCorrectPoints(pointField& p) const
{
    twoDPointCorrector::New(mesh_).correctPoints(p);
}


bool Foam::pointMeshMover::write() const
{
    return true;
}


// ************************************************************************* //
