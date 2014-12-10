/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "one.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cartesianCS, 0);
    addToRunTimeSelectionTable(coordinateSystem, cartesianCS, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cartesianCS::cartesianCS()
:
    coordinateSystem()
{}


Foam::cartesianCS::cartesianCS
(
    const coordinateSystem& cs
)
:
    coordinateSystem(cs)
{}


Foam::cartesianCS::cartesianCS
(
    const word& name,
    const coordinateSystem& cs
)
:
    coordinateSystem(name, cs)
{}


Foam::cartesianCS::cartesianCS
(
    const word& name,
    const point& origin,
    const coordinateRotation& cr
)
:
    coordinateSystem(name, origin, cr)
{}


Foam::cartesianCS::cartesianCS
(
    const word& name,
    const point& origin,
    const vector& axis,
    const vector& dirn
)
:
    coordinateSystem(name, origin, axis, dirn)
{}


Foam::cartesianCS::cartesianCS
(
    const word& name,
    const dictionary& dict
)
:
    coordinateSystem(name, dict)
{}


Foam::cartesianCS::cartesianCS
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    coordinateSystem(obr, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cartesianCS::~cartesianCS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::vector Foam::cartesianCS::localToGlobal
(
    const vector& local,
    bool translate
) const
{
    return coordinateSystem::localToGlobal(local, translate);
}


Foam::tmp<Foam::vectorField> Foam::cartesianCS::localToGlobal
(
    const vectorField& local,
    bool translate
) const
{
    return coordinateSystem::localToGlobal(local, translate);
}


Foam::vector Foam::cartesianCS::globalToLocal
(
    const vector& global,
    bool translate
) const
{
    return coordinateSystem::globalToLocal(global, translate);
}


Foam::tmp<Foam::vectorField> Foam::cartesianCS::globalToLocal
(
    const vectorField& global,
    bool translate
) const
{
    return coordinateSystem::globalToLocal(global, translate);
}


// ************************************************************************* //
