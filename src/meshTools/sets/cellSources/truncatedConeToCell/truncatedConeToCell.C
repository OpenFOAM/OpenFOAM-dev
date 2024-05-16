/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "truncatedConeToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(truncatedConeToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, truncatedConeToCell, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::truncatedConeToCell::combine(topoSet& set, const bool add) const
{
    const vector axis = point2_ - point1_;
    const scalar magAxisSquared = magSqr(axis);

    const pointField& ctrs = mesh_.cellCentres();

    forAll(ctrs, celli)
    {
        const vector d = ctrs[celli] - point1_;
        const scalar magD = d & axis;

        if ((magD > 0) && (magD < magAxisSquared))
        {
            const scalar d2 = (d & d) - sqr(magD) / magAxisSquared;
            const scalar radiusAtCell =
                radius1_ + (radius2_ - radius1_)*magD/(axis&axis);

            if (d2 <= sqr(radiusAtCell))
            {
                addOrDelete(set, celli, add);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::truncatedConeToCell::truncatedConeToCell
(
    const polyMesh& mesh,
    const vector& point1,
    const vector& point2,
    const scalar radius1,
    const scalar radius2
)
:
    topoSetSource(mesh),
    point1_(point1),
    point2_(point2),
    radius1_(radius1),
    radius2_(radius2)
{}


Foam::truncatedConeToCell::truncatedConeToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    point1_(dict.lookup<point>("point1", dimLength)),
    point2_(dict.lookup<point>("point2", dimLength)),
    radius1_(dict.lookup<scalar>("radius1", dimLength)),
    radius2_(dict.lookup<scalar>("radius2", dimLength))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::truncatedConeToCell::~truncatedConeToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::truncatedConeToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells with centre within cone, with point1 = "
            << point1_
            << ", point2 = " << point2_
            << " and radii = " << radius1_ << " to " << radius2_ << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells with centre within cone, with point1 = "
            << point1_
            << ", point2 = " << point2_
            << " and radii = " << radius1_ << " to " << radius2_ << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
