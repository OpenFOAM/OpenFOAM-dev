/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "cartesianCS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateSystems
{
    defineTypeNameAndDebug(cartesian, 0);
    addToRunTimeSelectionTable(coordinateSystem, cartesian, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::coordinateSystems::cartesian::localToGlobal
(
    const vector& local,
    bool translate
) const
{
    return coordinateSystem::localToGlobal(local, translate);
}


Foam::tmp<Foam::vectorField> Foam::coordinateSystems::cartesian::localToGlobal
(
    const vectorField& local,
    bool translate
) const
{
    return coordinateSystem::localToGlobal(local, translate);
}


Foam::vector Foam::coordinateSystems::cartesian::globalToLocal
(
    const vector& global,
    bool translate
) const
{
    return coordinateSystem::globalToLocal(global, translate);
}


Foam::tmp<Foam::vectorField> Foam::coordinateSystems::cartesian::globalToLocal
(
    const vectorField& global,
    bool translate
) const
{
    return coordinateSystem::globalToLocal(global, translate);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateSystems::cartesian::cartesian
(
    const word& name,
    const point& origin,
    const coordinateRotation& cr
)
:
    coordinateSystem(name, origin, cr)
{}


Foam::coordinateSystems::cartesian::cartesian
(
    const word& name,
    const point& origin,
    const vector& axis,
    const vector& dirn
)
:
    coordinateSystem(name, origin, axis, dirn)
{}


Foam::coordinateSystems::cartesian::cartesian
(
    const word& name,
    const dictionary& dict
)
:
    coordinateSystem(name, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordinateSystems::cartesian::~cartesian()
{}


// ************************************************************************* //
