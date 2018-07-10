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

#include "cylindricalCS.H"

#include "one.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylindricalCS::cylindricalCS(const bool inDegrees)
:
    coordinateSystem(),
    inDegrees_(inDegrees)
{}


Foam::cylindricalCS::cylindricalCS
(
    const coordinateSystem& cs,
    const bool inDegrees
)
:
    coordinateSystem(cs),
    inDegrees_(inDegrees)
{}


Foam::cylindricalCS::cylindricalCS
(
    const word& name,
    const coordinateSystem& cs,
    const bool inDegrees
)
:
    coordinateSystem(name, cs),
    inDegrees_(inDegrees)
{}


Foam::cylindricalCS::cylindricalCS
(
    const word& name,
    const point& origin,
    const coordinateRotation& cr,
    const bool inDegrees
)
:
    coordinateSystem(name, origin, cr),
    inDegrees_(inDegrees)
{}


Foam::cylindricalCS::cylindricalCS
(
    const word& name,
    const point& origin,
    const vector& axis,
    const vector& dirn,
    const bool inDegrees
)
:
    coordinateSystem(name, origin, axis, dirn),
    inDegrees_(inDegrees)
{}


Foam::cylindricalCS::cylindricalCS
(
    const word& name,
    const dictionary& dict
)
:
    coordinateSystem(name, dict),
    inDegrees_(dict.lookupOrDefault("degrees", true))
{}


Foam::cylindricalCS::cylindricalCS
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    coordinateSystem(obr, dict),
    inDegrees_(dict.lookupOrDefault("degrees", true))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cylindricalCS::~cylindricalCS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cylindricalCS::inDegrees() const
{
    return inDegrees_;
}


bool& Foam::cylindricalCS::inDegrees()
{
    return inDegrees_;
}


Foam::vector Foam::cylindricalCS::localToGlobal
(
    const vector& local,
    bool translate
) const
{
    scalar theta
    (
        local.y()*(inDegrees_ ? constant::mathematical::pi/180.0 : 1.0)
    );

    return coordinateSystem::localToGlobal
    (
        vector(local.x()*cos(theta), local.x()*sin(theta), local.z()),
        translate
    );
}


Foam::tmp<Foam::vectorField> Foam::cylindricalCS::localToGlobal
(
    const vectorField& local,
    bool translate
) const
{
    scalarField theta
    (
        local.component(vector::Y)
       *(inDegrees_ ? constant::mathematical::pi/180.0 : 1.0)
    );


    vectorField lc(local.size());
    lc.replace(vector::X, local.component(vector::X)*cos(theta));
    lc.replace(vector::Y, local.component(vector::X)*sin(theta));
    lc.replace(vector::Z, local.component(vector::Z));

    return coordinateSystem::localToGlobal(lc, translate);
}


Foam::vector Foam::cylindricalCS::globalToLocal
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
        atan2
        (
            lc.y(),
            lc.x()
        )*(inDegrees_ ? 180.0/constant::mathematical::pi : 1.0),
        lc.z()
    );
}


Foam::tmp<Foam::vectorField> Foam::cylindricalCS::globalToLocal
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
        atan2
        (
            lc.component(vector::Y),
            lc.component(vector::X)
        )*(inDegrees_ ? 180.0/constant::mathematical::pi : 1.0)
    );

    result.replace(vector::Z, lc.component(vector::Z));

    return tresult;
}


// ************************************************************************* //
