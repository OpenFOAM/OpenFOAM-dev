/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "patchDistanceToCell.H"
#include "patchWave.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchDistanceToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, patchDistanceToCell, word);
    addToRunTimeSelectionTable(topoSetSource, patchDistanceToCell, istream);
}


Foam::topoSetSource::addToUsageTable Foam::patchDistanceToCell::usage_
(
    patchDistanceToCell::typeName,
    "\n    Usage: patchDistanceToCell (<patches>) distance\n\n"
    "    Select cells that are below a distance from a list of patches\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchDistanceToCell::combine(topoSet& set, const bool add) const
{
    const patchWave pw
    (
        mesh_,
        mesh_.boundaryMesh().patchSet(patches_)
    );

    forAll(pw.distance(), celli)
    {
        if (pw.distance()[celli] < distance_)
        {
            addOrDelete(set, celli, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDistanceToCell::patchDistanceToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    patches_
    (
        dict.found("patches")
      ? dict.lookup<wordReList>("patches")
      : wordReList(1, dict.lookup<wordRe>("patch"))
    ),
    distance_(dict.lookup<scalar>("distance"))
{}


Foam::patchDistanceToCell::patchDistanceToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    patches_(),
    distance_()
{
    token firstToken(is);
    is.putBack(firstToken);

    if (firstToken == token::BEGIN_LIST)
    {
        is >> patches_;
    }
    else
    {
        patches_ = wordReList(1, word(is));
    }

    is >> distance_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchDistanceToCell::~patchDistanceToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchDistanceToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells a distance less than " << distance_
            << " from patches " << patches_ << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells a distance less than " << distance_
            << " from patches " << patches_ << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
