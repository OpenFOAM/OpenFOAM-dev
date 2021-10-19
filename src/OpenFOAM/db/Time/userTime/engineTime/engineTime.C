/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "engineTime.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace userTimes
{
    defineTypeNameAndDebug(engine, 0);
    addToRunTimeSelectionTable(userTime, engine, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::userTimes::engine::engine(const dictionary& controlDict)
:
    userTime(controlDict),
    rpm_("rpm", dimless/dimTime, 0)
{
    read(controlDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::userTimes::engine::~engine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dimensionedScalar& Foam::userTimes::engine::rpm() const
{
    return rpm_;
}


Foam::scalar Foam::userTimes::engine::userTimeToTime
(
    const scalar theta
) const
{
    return theta/(6*rpm_.value());
}


Foam::scalar Foam::userTimes::engine::timeToUserTime(const scalar t) const
{
    return t*(6*rpm_.value());
}


Foam::word Foam::userTimes::engine::unit() const
{
    return "CAD";
}


bool Foam::userTimes::engine::read(const dictionary& controlDict)
{
    rpm_.read(dict(controlDict));
    return true;
}


// ************************************************************************* //
