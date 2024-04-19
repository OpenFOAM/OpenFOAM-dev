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

#include "cylindricalCS.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateSystems
{
    defineTypeNameAndDebug(cylindrical, 0);
    addToRunTimeSelectionTable(coordinateSystem, cylindrical, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::coordinateSystems::cylindrical::localToGlobal
(
    const vector& local,
    bool translate
) const
{
    const scalar theta = local.y();

    return coordinateSystem::localToGlobal
    (
        vector(local.x()*cos(theta), local.x()*sin(theta), local.z()),
        translate
    );
}


Foam::tmp<Foam::vectorField> Foam::coordinateSystems::cylindrical::localToGlobal
(
    const vectorField& local,
    bool translate
) const
{
    const scalarField theta(local.component(vector::Y));

    vectorField lc(local.size());
    lc.replace(vector::X, local.component(vector::X)*cos(theta));
    lc.replace(vector::Y, local.component(vector::X)*sin(theta));
    lc.replace(vector::Z, local.component(vector::Z));

    return coordinateSystem::localToGlobal(lc, translate);
}


Foam::vector Foam::coordinateSystems::cylindrical::globalToLocal
(
    const vector& global,
    bool translate
) const
{
    const vector lc
    (
        coordinateSystem::globalToLocal(global, translate)
    );

    return vector
    (
        sqrt(sqr(lc.x()) + sqr(lc.y())),
        atan2(lc.y(), lc.x()),
        lc.z()
    );
}


Foam::tmp<Foam::vectorField> Foam::coordinateSystems::cylindrical::globalToLocal
(
    const vectorField& global,
    bool translate
) const
{
    const vectorField lc
    (
        coordinateSystem::globalToLocal(global, translate)
    );

    tmp<vectorField> tresult(new vectorField(lc.size()));
    vectorField& result = tresult.ref();

    result.replace
    (
        vector::X,
        sqrt(sqr(lc.component(vector::X)) + sqr(lc.component(vector::Y)))
    );

    result.replace
    (
        vector::Y,
        atan2(lc.component(vector::Y), lc.component(vector::X))
    );

    result.replace(vector::Z, lc.component(vector::Z));

    return tresult;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateSystems::cylindrical::cylindrical
(
    const word& name,
    const point& origin,
    const coordinateRotation& cr
)
:
    coordinateSystem(name, origin, cr)
{}


Foam::coordinateSystems::cylindrical::cylindrical
(
    const word& name,
    const point& origin,
    const vector& axis,
    const vector& dirn
)
:
    coordinateSystem(name, origin, axis, dirn)
{}


Foam::coordinateSystems::cylindrical::cylindrical
(
    const word& name,
    const dictionary& dict
)
:
    coordinateSystem(name, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordinateSystems::cylindrical::~cylindrical()
{}


// ************************************************************************* //
