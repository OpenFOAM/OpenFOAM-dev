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

#include "cylinderAnnulusToFace.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cylinderAnnulusToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, cylinderAnnulusToFace, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cylinderAnnulusToFace::combine(topoSet& set, const bool add) const
{
    const vector axis = point2_ - point1_;
    const scalar orad2 = sqr(outerRadius_);
    const scalar irad2 = sqr(innerRadius_);
    const scalar magAxis2 = magSqr(axis);

    const pointField& ctrs = mesh_.faceCentres();

    forAll(ctrs, facei)
    {
        vector d = ctrs[facei] - point1_;
        scalar magD = d & axis;

        if ((magD > 0) && (magD < magAxis2))
        {
            scalar d2 = (d & d) - sqr(magD)/magAxis2;
            if ((d2 < orad2) && (d2 > irad2))
            {
                addOrDelete(set, facei, add);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylinderAnnulusToFace::cylinderAnnulusToFace
(
    const polyMesh& mesh,
    const vector& point1,
    const vector& point2,
    const scalar outerRadius,
    const scalar innerRadius
)
:
    topoSetSource(mesh),
    point1_(point1),
    point2_(point2),
    outerRadius_(outerRadius),
    innerRadius_(innerRadius)
{}


Foam::cylinderAnnulusToFace::cylinderAnnulusToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    point1_(dict.lookupBackwardsCompatible<point>({"point1", "p1"})),
    point2_(dict.lookupBackwardsCompatible<point>({"point2", "p2"})),
    outerRadius_(dict.lookup<scalar>("outerRadius")),
    innerRadius_(dict.lookup<scalar>("innerRadius"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cylinderAnnulusToFace::~cylinderAnnulusToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cylinderAnnulusToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding faces with centre within cylinder annulus,"
            << " with point1 = "
            << point1_ << ", point2 = " << point2_ << " and outer radius = "
            << outerRadius_ << " and inner radius = " << innerRadius_ << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing faces with centre within cylinder, with point1 = "
            << point1_ << ", point2 = " << point2_ << " and outer radius = "
            << outerRadius_ << " and inner radius " << innerRadius_ << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
