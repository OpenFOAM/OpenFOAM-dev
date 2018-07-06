/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "cylinderAnnulusToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cylinderAnnulusToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, cylinderAnnulusToCell, word);
    addToRunTimeSelectionTable(topoSetSource, cylinderAnnulusToCell, istream);
}


Foam::topoSetSource::addToUsageTable Foam::cylinderAnnulusToCell::usage_
(
    cylinderAnnulusToCell::typeName,
    "\n    Usage: cylinderAnnulusToCell (p1X p1Y p1Z) (p2X p2Y p2Z)"
    " outerRadius innerRadius\n\n"
    "    Select all cells with cell centre within bounding cylinder annulus\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cylinderAnnulusToCell::combine(topoSet& set, const bool add) const
{
    const vector axis = p2_ - p1_;
    const scalar orad2 = sqr(outerRadius_);
    const scalar irad2 = sqr(innerRadius_);
    const scalar magAxis2 = magSqr(axis);

    const pointField& ctrs = mesh_.cellCentres();

    forAll(ctrs, celli)
    {
        vector d = ctrs[celli] - p1_;
        scalar magD = d & axis;

        if ((magD > 0) && (magD < magAxis2))
        {
            scalar d2 = (d & d) - sqr(magD)/magAxis2;
            if ((d2 < orad2) && (d2 > irad2))
            {
                addOrDelete(set, celli, add);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylinderAnnulusToCell::cylinderAnnulusToCell
(
    const polyMesh& mesh,
    const vector& p1,
    const vector& p2,
    const scalar outerRadius,
    const scalar innerRadius
)
:
    topoSetSource(mesh),
    p1_(p1),
    p2_(p2),
    outerRadius_(outerRadius),
    innerRadius_(innerRadius)
{}


Foam::cylinderAnnulusToCell::cylinderAnnulusToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    p1_(dict.lookup("p1")),
    p2_(dict.lookup("p2")),
    outerRadius_(readScalar(dict.lookup("outerRadius"))),
    innerRadius_(readScalar(dict.lookup("innerRadius")))
{}


Foam::cylinderAnnulusToCell::cylinderAnnulusToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    p1_(checkIs(is)),
    p2_(checkIs(is)),
    outerRadius_(readScalar(checkIs(is))),
    innerRadius_(readScalar(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cylinderAnnulusToCell::~cylinderAnnulusToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cylinderAnnulusToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells with centre within cylinder annulus,"
            << " with p1 = "
            << p1_ << ", p2 = " << p2_ << " and outer radius = " << outerRadius_
        << " and inner radius = " << innerRadius_
        << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells with centre within cylinder, with p1 = "
            << p1_ << ", p2 = " << p2_ << " and outer radius = " << outerRadius_
        << " and inner radius " << innerRadius_
        << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
