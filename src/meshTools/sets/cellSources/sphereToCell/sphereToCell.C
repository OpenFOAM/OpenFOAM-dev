/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "sphereToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sphereToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, sphereToCell, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sphereToCell::combine(topoSet& set, const bool add) const
{
    const pointField& ctrs = mesh_.cellCentres();

    const scalar radSquared = radius_*radius_;

    forAll(ctrs, celli)
    {
        scalar offset = magSqr(centre_ - ctrs[celli]);
        if (offset <= radSquared)
        {
            addOrDelete(set, celli, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sphereToCell::sphereToCell
(
    const polyMesh& mesh,
    const vector& centre,
    const scalar radius
)
:
    topoSetSource(mesh),
    centre_(centre),
    radius_(radius)
{}


Foam::sphereToCell::sphereToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    centre_(dict.lookup<point>("centre", dimLength)),
    radius_(dict.lookup<scalar>("radius", dimLength))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sphereToCell::~sphereToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sphereToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells with centre within sphere, with centre = "
            << centre_ << " and radius = " << radius_ << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells with centre within sphere, with centre = "
            << centre_ << " and radius = " << radius_ << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
