/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    addToRunTimeSelectionTable(topoSetSource, sphereToCell, istream);
}


Foam::topoSetSource::addToUsageTable Foam::sphereToCell::usage_
(
    sphereToCell::typeName,
    "\n    Usage: sphereToCell (centreX centreY centreZ) radius\n\n"
    "    Select all cells with cellCentre within bounding sphere\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sphereToCell::combine(topoSet& set, const bool add) const
{
    const pointField& ctrs = mesh_.cellCentres();

    const scalar radSquared = radius_*radius_;

    forAll(ctrs, cellI)
    {
        scalar offset = magSqr(centre_ - ctrs[cellI]);
        if (offset <= radSquared)
        {
            addOrDelete(set, cellI, add);
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
    centre_(dict.lookup("centre")),
    radius_(readScalar(dict.lookup("radius")))
{}


// Construct from Istream
Foam::sphereToCell::sphereToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    centre_(checkIs(is)),
    radius_(readScalar(checkIs(is)))
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
