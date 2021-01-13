/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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

#include "cylinderToFace.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cylinderToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, cylinderToFace, word);
    addToRunTimeSelectionTable(topoSetSource, cylinderToFace, istream);
}


Foam::topoSetSource::addToUsageTable Foam::cylinderToFace::usage_
(
    cylinderToFace::typeName,
    "\n    Usage: cylinderToFace (point1X point1Y point1Z)"
    " (point2X point2Y point2Z) radius\n\n"
    "    Select all faces with face centre within bounding cylinder\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cylinderToFace::combine(topoSet& set, const bool add) const
{
    const vector axis = point2_ - point1_;
    const scalar rad2 = sqr(radius_);
    const scalar magAxis2 = magSqr(axis);

    const pointField& ctrs = mesh_.faceCentres();

    forAll(ctrs, facei)
    {
        vector d = ctrs[facei] - point1_;
        scalar magD = d & axis;

        if ((magD > 0) && (magD < magAxis2))
        {
            scalar d2 = (d & d) - sqr(magD)/magAxis2;
            if (d2 < rad2)
            {
                addOrDelete(set, facei, add);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylinderToFace::cylinderToFace
(
    const polyMesh& mesh,
    const vector& point1,
    const vector& point2,
    const scalar radius
)
:
    topoSetSource(mesh),
    point1_(point1),
    point2_(point2),
    radius_(radius)
{}


Foam::cylinderToFace::cylinderToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    point1_(dict.lookupBackwardsCompatible<point>({"point1", "p1"})),
    point2_(dict.lookupBackwardsCompatible<point>({"point2", "p2"})),
    radius_(dict.lookup<scalar>("radius"))
{}


Foam::cylinderToFace::cylinderToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    point1_(checkIs(is)),
    point2_(checkIs(is)),
    radius_(readScalar(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cylinderToFace::~cylinderToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cylinderToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding faces with centre within cylinder, with point1 = "
            << point1_ << ", point2 = " << point2_ << " and radius = "
            << radius_ << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing faces with centre within cylinder, with point1 = "
            << point1_ << ", point2 = " << point2_ << " and radius = "
            << radius_ << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
